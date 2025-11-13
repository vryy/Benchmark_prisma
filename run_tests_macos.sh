#!/bin/zsh
reset
# # recompile the code
# cd $HOME/kratos_builds/cmake_build_release
# sh do_make.sh
# cd -
# ignore system error warnings
ulimit -c 0
# run tests
export PATH="$HOME/opt/cmake-3.30.5/bin:$HOME/opt/openmpi-3.1.2/bin:/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin:/opt/pmk/env/global/bin:/opt/X11/bin:/Library/Apple/usr/bin:/Applications/iTerm.app/Contents/Resources/utilities"
export BENCHMARK_PRISMA=$PWD
export KRATOS_ROOT_PATH=$HOME/kratos_bundle/kratos_bcn2
export KRATOS_INSTALL_PREFIX=$HOME/kratos_bundle/kratos_libs/Release314
export PYTHONPATH=$KRATOS_INSTALL_PREFIX:$KRATOS_INSTALL_PREFIX/libs
export LD_LIBRARY_PATH=$KRATOS_INSTALL_PREFIX/libs:$HOME/opt/boost_1_86_0/lib:$HOME/kratos_bundle/kratos_libs/Release314/libs:$HOME/opt/hdf5-1.8.21/lib:$HOME/opt/qhull-2020.2/lib:$HOME/opt/adol-c-2.7.2/lib:$HOME/opt/mmg-5.5.2/lib:$HOME/opt/hdf5-1.8.21/lib:$HOME/opt/p4est-vryy/lib:$HOME/opt/petsc-3.22.5/lib:$HOME/opt/petsc-3.22.5/lib:$HOME/opt/openmpi-3.1.2/lib:$HOME/opt/cgal-6.1/lib:$HOME/opt/gmp-6.3.0/lib:$HOME/opt/T-SPLINE/lib/x86/release:$HOME/opt/cddlib-094j/lib:$HOME/opt/mpfr-3.1.2/lib
export OMP_NUM_THREADS=1
export PY_COMMAND="python"
mkdir -p ./ztest_logs
if [ "$#" -eq 0 ] ; then
    echo "Run all tests"
    export OUTPUT="./ztest_logs/$(whoami)_$(hostname)_"$(date +%Y_%m_%d)".log"
    python run_tests.py $PY_COMMAND | tee $OUTPUT
elif [ "$#" -ge 1 ] ; then
    echo "Run tests with arguments:" $@
    export OUTPUT="./ztest_logs/$(whoami)_$(hostname)_"$(date +%Y_%m_%d)"-$1.log"
    python run_tests.py $PY_COMMAND $@ | tee $OUTPUT
fi
