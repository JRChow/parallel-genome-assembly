set -e
export UPCXX_SEGMENT_MB=4000
export GASNET_MAX_SEGSIZE=4G
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC -DKMER_LEN=51 .. && cmake --build .
DATA_DIR=$SCRATCH/my_datasets
n_node=4
n_proc=272
salloc -N $n_node -A mp309 -t 30:00 -q debug --qos=interactive -C knl srun -N $n_node -n $n_proc ./kmer_hash $DATA_DIR/human-chr14-synthetic.txt test
echo "Checking correctness..."
cat test*.dat | sort > my_solution.txt
diff my_solution.txt $DATA_DIR/human-chr14-synthetic_solution.txt

