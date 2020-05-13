#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include "butil.hpp"

struct HashMap {
    // std::vector<int> used;
    size_t global_hashmap_size;  // Total number of entries in global hash map
    size_t offset;   // Number of entries between every two processes
    size_t my_base;  // Number of entries before me belonging to other processes
    size_t my_size;  // Number of entries belonging to me

    int num_proc;  // Total number of processes
    int my_rank;

    /* Suggested implementation on class website: vector of global_ptrs that
     * point to arrays allocated in the shared segment on each rank. */
    std::vector <upcxx::global_ptr<kmer_pair>> global_map_ptrs;
    std::vector <upcxx::global_ptr<int>> is_occupied_slots;
    // Atomic domain used to access the above occupancy list
    upcxx::atomic_domain<int> dom = upcxx::atomic_domain<int>({upcxx::atomic_op::compare_exchange,
                                                               upcxx::atomic_op::load});

    size_t size() const noexcept;  // Getter for my_size

    HashMap(size_t global_hashmap_size);  // Constructor
    ~HashMap();                           // Destructor

    // Most important functions: insert and retrieve k-mers from the hash table.
    bool insert(const kmer_pair &kmer);

    bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

    // Helper functions
    inline uint64_t get_global_idx(uint64_t hashcode);

    inline upcxx::global_ptr <kmer_pair> get_ptr_to_global_elem(uint64_t global_idx);

    inline upcxx::global_ptr<int> get_ptr_to_global_occupancy(uint64_t global_idx);

    inline int CAS_on_global_occupancy(uint64_t global_idx);

    inline bool is_slot_empty(uint64_t global_idx);

    inline void next_slot(uint64_t* global_idx_ptr);
};

// [Jingran] Constructor
HashMap::HashMap(size_t global_hashmap_size) {
    // Initialize member variables
    this->global_hashmap_size = global_hashmap_size;
    num_proc = upcxx::rank_n();
    my_rank = upcxx::rank_me();
    // Equivalent to ceil(global_hashmap_size / num_proc)
    offset = (global_hashmap_size + num_proc - 1) / num_proc;
    my_base = offset * my_rank;
    my_size = offset;
    // Corner case: last process might have a different size!
    if (my_rank == num_proc - 1) {
        my_size = global_hashmap_size - (num_proc - 1) * offset;
    }
    // Prepare vector of global pointers
    global_map_ptrs = std::vector < upcxx::global_ptr < kmer_pair >> (num_proc);
    is_occupied_slots = std::vector < upcxx::global_ptr < int >> (num_proc);

    // FIXME: do we really need a barrier here?
    upcxx::barrier();

    // Each process builds its part of the global hashmap and broadcast to all
    upcxx::global_ptr <kmer_pair> my_part = upcxx::new_array<kmer_pair>(my_size);
    upcxx::global_ptr<int> my_occupancy = upcxx::new_array<int>(my_size);
    // Broadcast my part to all and receive other parts from all
    for (int p = 0; p < num_proc; ++p) {
        global_map_ptrs[p] = upcxx::broadcast(my_part, p).wait();
        is_occupied_slots[p] = upcxx::broadcast(my_occupancy, p).wait();
    }

    // FIXME: do we really need a barrier here?
    upcxx::barrier();
}

// [Jingran] Destructor
HashMap::~HashMap() { this->dom.destroy(); }

// [Daniel] Insert k-mer pair into the hash map using linear probing
bool HashMap::insert(const kmer_pair &kmer) {
    uint64_t g_idx = get_global_idx(kmer.hash());
    uint64_t initial_g_idx = g_idx;
    // Linear probing
    while (CAS_on_global_occupancy(g_idx) != 0) {
        next_slot(&g_idx);
        if (initial_g_idx == g_idx)
            return false;
    }
    upcxx::rput(kmer, get_ptr_to_global_elem(g_idx)).wait();
    return true;
}

// [Jingran] Finding item in global hashmap with key_kmer and store pair in val_kmer
bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
    uint64_t g_idx = get_global_idx(key_kmer.hash());
    uint64_t initial_g_idx = g_idx;
    do {
        // If encountering an empty slot, we've failed.
        if (is_slot_empty(g_idx)) break;
        // Retrieve element from global hash map
        kmer_pair item = upcxx::rget(get_ptr_to_global_elem(g_idx)).wait();
        // If keys match, we've found our target
        if (item.kmer == key_kmer) {
            val_kmer = item;
            return true;
        }
        next_slot(&g_idx);
    } while (g_idx != initial_g_idx);
    // Key not found
    return false;
}

// Given a hashcode, determine the global index in the hash map
inline uint64_t HashMap::get_global_idx(uint64_t hashcode) { return hashcode % global_hashmap_size; }

// Given a global index, get a global pointer to the element it refers to
inline upcxx::global_ptr <kmer_pair> HashMap::get_ptr_to_global_elem(uint64_t global_idx) {
    return global_map_ptrs[global_idx / offset] + (global_idx % offset);
}

// Given a global index, get a global pointer to the occupancy it refers to
inline upcxx::global_ptr<int> HashMap::get_ptr_to_global_occupancy(uint64_t global_idx) {
    return is_occupied_slots[global_idx / offset] + (global_idx % offset);
}

// Compare-and-Exchange on global occupancy array
inline int HashMap::CAS_on_global_occupancy(uint64_t global_idx) {
    /* If this slot is available (equals 0), occupy it (make it 1) and return 0,
     * indicating we've found an empty place. If this slot is occupied (equals 1),
     * do nothing and return 1, indicating it's taken. */
    return dom.compare_exchange(get_ptr_to_global_occupancy(global_idx),
                                0, 1, std::memory_order_relaxed).wait();
}

// Check if the slot pointed to is non-empty
inline bool HashMap::is_slot_empty(uint64_t global_idx) {
    return dom.load(get_ptr_to_global_occupancy(global_idx), std::memory_order_relaxed).wait() == 0;
}

// Go to the next slot in the hash map (modulo wrapped)
inline void HashMap::next_slot(uint64_t* global_idx_ptr) {
    *global_idx_ptr += 1;
    *global_idx_ptr %= global_hashmap_size;  // Wrap around
}
