#!/bin/bash -x

narr=400
nt=100
tsc=2.593906
np=64
host=node923

mpicc -O2 -Wall -o jacobi2d5p_cgt_tf.x ./jacobi2d5p_tf.c  -DTIMING -DUSE_CGT -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas
mpicc -O2 -Wall -o jacobi2d5p_papi_tf.x ./jacobi2d5p_tf.c  -DTIMING -DUSE_PAPI -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi                               
mpicc -O2 -Wall -o jacobi2d5p_papix6_tf.x ./jacobi2d5p_tf.c  -DTIMING -DUSE_PAPIX6 -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi                           
mpicc -O2 -Wall -o jacobi2d5p_likwid_tf.x ./jacobi2d5p_tf.c  -DTIMING -DUSE_LIKWID -DLIKWID_PERFMON -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -llikwid

mpicc -O2 -Wall -o jacobi2d5p_cgt.x ./jacobi2d5p.c  -DTIMING -DUSE_CGT -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas
mpicc -O2 -Wall -o jacobi2d5p_papi.x ./jacobi2d5p.c  -DTIMING -DUSE_PAPI -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi
mpicc -O2 -Wall -o jacobi2d5p_papix6.x ./jacobi2d5p.c  -DTIMING -DUSE_PAPIX6 -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi                        
mpicc -O2 -Wall -o jacobi2d5p_likwid.x ./jacobi2d5p.c  -DTIMING -DUSE_LIKWID -DLIKWID_PERFMON -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -llikwid  

for m in cgt papi papix6
do
    mpirun --map-by core --bind-to core -np ${np} ./jacobi2d5p_${m}.x $narr $tsc
    rm  -r jacobi2d5p${narr}n${nt}t_${m}_${host}
    mkdir jacobi2d5p${narr}n${nt}t_${m}_${host}
    mv ./*.csv jacobi2d5p${narr}n${nt}t_${m}_${host}
done

likwid-mpirun -np ${np} -g L3 -m ./jacobi2d5p_likwid.x $narr $tsc
rm  -r jacobi2d5p${narr}n${nt}t_likwid_${host}
mkdir jacobi2d5p${narr}n${nt}t_likwid_${host}
mv ./*.csv jacobi2d5p${narr}n${nt}t_likwid_${host}

