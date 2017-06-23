#!/bin/bash -l
#SBATCH --job-name="test"
#SBATCH -p regular
#SBATCH --nodes=1
#SBATCH -t 00:40:00
#SBATCH -o test_nwchem.o%j

export OMP_NUM_THREADS=48
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export TURBODIR=$HOME/apps/lrsh PARA_ARCH=SMP
source $TURBODIR/Config_turbo_env
env> env.out

actual -r
kdg dft
sed -i '1 i $scftol 1d-12\n$profile\n$dft\nfunctional pbe\nreference' control
dscf_smp > dscf.pbe.gen
kdg scfiterlimit
sed -i '1 i $scfiterlimit 1' control
cpc ini

for i in pbe pbe0 tpss tpssh b3-lyp b-p; do
    cp ini/* .
    kdg dft
    sed -i '1 i $dft\nfunctional '$i'\nreference' control
    dscf_smp > dscf.$i
done

cp ini/* .
kdg dft
sed -i '1 i $rshexp 2\n$nlbeta 0d0\n$dft\nfunctional lrsh\ngridsize 5\n$nonlscreen' control
dscf_smp> dscf.lrsh

echo Job Done
