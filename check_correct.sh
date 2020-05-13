set -e
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC -DKMER_LEN=19 .. && cmake --build .
DATA_DIR=$SCRATCH/my_datasets
# DATA_DIR=$SCRATCH/my_datasets/smaller
n_proc=68
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C knl srun -N 1 -n $n_proc ./kmer_hash $DATA_DIR/test.txt test
echo "Checking correctness..."
cat test*.dat | sort > my_solution.txt
diff my_solution.txt $DATA_DIR/test_solution.txt

