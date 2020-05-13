set -e
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC -DKMER_LEN=19 .. && cmake --build .
DATA_DIR=$SCRATCH/my_datasets
# DATA_DIR=$SCRATCH/my_datasets/smaller
N_NODE=1
N_PROC=`expr $N_NODE \* 68`
MODE=verbose
salloc -N $N_NODE -A mp309 -t 10:00 -q debug --qos=interactive -C knl srun -N $N_NODE -n $N_PROD ./kmer_hash $DATA_DIR/test.txt $MODE
if [ $MODE == "test" ]; then
    echo "Checking correctness..."
    cat test*.dat | sort > my_solution.txt
    diff my_solution.txt $DATA_DIR/test_solution.txt
fi
