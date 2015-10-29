#! /bin/bash
# execute from TESTSET root: 
#    bash unittest/test-all-the-thingies.sh [method] [basis] [testset]
# method, basis, testset are required, 
## USE LRSH FOR INITIAL TESTING HAS CORRECT DEFINE INPUT
# method defines which file will be used for setup i.e. setupfiles/setup.method 
# runs smp version of dscf just for speed sake, though sequential is probably safer
# testset variable can either be a single test or an array of tests:
#    bash unittest/test-all-the-thingies.sh [method] [basis] (g21ip g21ea)
# if testset variable is set to "all" then it will parse and set all testset dirs 
meth=$1
basis=$2
testset=$3

if [ "$testset" == "all" ]; then
    testset=$(ls | grep -v setupfiles | grep -v README | grep -v unittest)
fi

export TURBODIR=/modfac/apps/TURBOMOLE_6.5_modfac PARA_ARCH=SMP 
export OMP_NUM_THREADS=$(cat /proc/cpuinfo | grep proc | wc -l)
source $TURBODIR/Config_turbo_env

for i in $testset; do
    ./setupfiles/setupTESTSET $meth $basis $i
    cd $i/$basis.$meth
    for j in $(ls); do
        cd $j
        dscf_smp > dscf.hf
        sed -i "1 i \$dft\nfunctional pbe\ngridsize m3" control
        dscf_smp > dscf.pbe
        sed -i "s/pbe/tpss/" control
        dscf_smp > dscf.tpss
        sed -i "s/tpss/b3\-lyp/" control
        dscf_smp > dscf.b3-lyp
        cd ..
    done
    python ../analyze*.py > res.out
    cd ../..
done
