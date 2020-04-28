#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    // std::vector<int> used;
    size_t global_hashmap_size;  // Total number of entries in global hash map
    size_t offset;   // Number of entries between every two processes
    size_t my_base;  // Number of entires before me belonging to other processes
    size_t my_size;  // Number of entries belonging to me

    int num_proc;  // Total number of processes
    int my_rank;

    /* Suggested implementation on class website: vector of global_ptrs that
     * point to arrays allocated in the shared segment on each rank. */
    std::vector<upcxx::global_ptr<kmer_pair>> global_map_ptrs;
    std::vector<upcxx::global_ptr<int>> is_occupied_slots;

    size_t size() const noexcept;  // Getter for my_size

    HashMap(size_t global_hashmap_size);  // Constructor

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // TODO: helper functions
    uint64_t hash_to_global_slot(uint64_t hash);
    bool check_slot_available(uint64_t slot);


    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

// [Jingran] Constructor
HashMap::HashMap(size_t global_hashmap_size) {
  // Initialize member variables
  this.global_hashmap_size = global_hashmap_size;
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

  global_map_ptrs = std::vector<upcxx::global_ptr<kmer_pair>>(num_proc);
  is_occupied_slots = std::vector<upcxx::global_ptr<int>>(num_proc);

  // FIXME: do we really need a barrier here?
  upcxx::barrier();

  // Each process builds its part of the global hashmap and broadcast to all
  upcxx::global_ptr<kmer_pair> my_part = new_array<kmer_pair>(my_size);
  upcxx::global_ptr<int> my_occupancy = new_array<int>(my_size);
  // Broadcast my part to all and receive other parts from all
  for (int p = 0; p < num_proc; ++p) {
    global_map_ptrs[p] = upcxx::broadcast(my_part, my_rank).wait();
    is_occupied_slots[p] = upcxx::broadcast(my_occupancy, my_rank).wait();
  }

  // FIXME: do we really need a barrier here?
  upcxx::barrier();
}

/**
bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        uint64_t slot = (hash + probe++) % size();
        success = request_slot(slot);
        if (success) {
            write_slot(slot, kmer);
        }
    } while (!success && probe < size());
    return success;
}
**/

// [Daniel]
bool HashMap::insert(const kmer_pair& kmer) {  // TODO
  // Compute the global slot
  unit64_t g_slot = kmer.hash() % global_hashmap_size;
  uint64_t initial_g_slot = g_slot;

  upcxx::global_ptr slot_ptr = is_occupied_slots[g_slot / offset] + (g_slot % offset);

  while(atomic_op::compare_exchange(slot_ptr, 0, 1) != 0) {
    g_slot = (g_slot+1) % global_hashmap_size;
    slot_ptr = is_occupied_slots[g_slot / offset] + (g_slot % offset);

    if (initial_g_slot == g_slot) {
      return false;
    }
  }

  upcxx::rput(kmer, slot_ptr).wait();
  return true;
}

// [Jingran] Finding item in global hashmap with key_kmer and store pair in val_kmer
bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {  // TODO
  uint64_t g_slot = key_kmer.hash() % global_hashmap_size;  // Global slot
  uint64_t initial_g_slot = g_slot;
  do {
    // Retrieve element from global hashmap
    int which_proc = g_slot / this.offset;
    uint64_t l_slot = g_slot % this.offset;  // Local slot
    // FIXME: probably don't need lock if we separate inserts from finds?
    kmer_pair item = upcxx::rget(global_map_ptrs[which_proc] + l_slot).wait();
    // If keys match, we've found our target
    if (item.kmer == key_kmer) {
      val_kmer = item;
      return true;
    }
    g_slot = (g_slot + 1) % global_hashmap_size;  // Linear probing
  } while (g_slot != initial_g_slot);

  return false;
}

// [Daniel]
/**
bool HashaMap::check_slot_available(uint64_t slot){
  // FIXME: Why not is_occupied_slots[slot] != 0?
  if (is_occupied_slots[slot]==0){
    return False;
  } else{
    return True;
  }
}
**/
//how to add lock, on that particular memory or only one entry?
//Each atomic operation works on a global pointer.
// atomic_domain<int64_t> dom({atomic_op::load, atomic_op::min, atomic_op::fetch_add, atomic_op::compare_exchange});
// uint_64 HashMap::find_open_slot(uint_64 slot_pos){
//   uint_64 temp_pos = slot_pos;
//   upcxx::global_ptr slot_ptr;
//   do{
//     temp_pos++;
//     slot_ptr = &is_occupied_slots[0]+temp_pos;
//   }while(atomic_op::compare_exchange(slot_ptr, 0, 1)!=0 && temp_pos<mysize) //mysize or offset
//
//   return temo_pos;
// }

bool HashMap::slot_used(uint64_t slot) {
  // return used[slot] != 0;
}

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) {
  // data[slot] = kmer;
}

kmer_pair HashMap::read_slot(uint64_t slot) {
  // return data[slot];
}

bool HashMap::request_slot(uint64_t slot) {
    // if (used[slot] != 0) {
    //     return false;
    // } else {
    //     used[slot] = 1;
    //     return true;
    // }
}

size_t HashMap::size() const noexcept { // Jingran: I don't see why we should keep this getter lol
  return my_size;
}

// [Jingran] Convert hash value to slot position in global hashmap
// uint64_t hash_to_global_slot(uint64_t hash) {
//   int which_proc = hash / this.offset;
//   uint64_t l_slot = hash % this.offset;  // Local slot
// }
