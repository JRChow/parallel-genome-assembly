set -e
cmake --build .
DATA_DIR='/project/projectdirs/mp309/cs267-spr2020/hw3-datasets'
n_proc=60
# salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C knl srun -N 1 -n $n_proc ./kmer_hash $DATA_DIR/smaller/tiny.txt test
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C knl srun -N 1 -n $n_proc ./kmer_hash $DATA_DIR/test.txt test
# echo "Checking correctness..."
cat test*.dat | sort > my_solution.txt
# diff my_solution.txt $DATA_DIR/smaller/tiny_solution.txt
diff my_solution.txt $DATA_DIR/test_solution.txt

