set -e
export UPCXX_SEGMENT_MB=4000
export GASNET_MAX_SEGSIZE=4G
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC -DKMER_LEN=51 .. && cmake --build .
DATA_DIR=$SCRATCH/my_datasets
N_NODE=10
N_PROC=`expr $N_NODE \* 68`
MODE=verbose
salloc -N $N_NODE -A mp309 -t 10:00 -q debug --qos=interactive -C knl srun -N $N_NODE -n $N_PROC ./kmer_hash $DATA_DIR/human-chr14-synthetic.txt $MODE
if [ $MODE == "test" ]; then
    echo "Checking correctness..."
    cat test*.dat | sort > my_solution.txt
    diff my_solution.txt $DATA_DIR/human-chr14-synthetic_solution.txt
fi
