
rm -rf a.m b.m
python random_float_matrix.py 1024 1024 >> a.m
python random_float_matrix.py 1024 1024 >> b.m

make mpi
make sequential

mpiexec -n 32 ./bin/mpi a.m b.m
./bin/seq a.m b.m