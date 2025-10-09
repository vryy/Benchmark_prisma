#!/bin/bash -l
reset

# define number of threads used for testing
NUM_THREADS=1

# test host folder
HOST_DIR="${HOME}/kratos_bundle/Benchmark_kratos"

# test run folder
SIMULATION_DIR="python3"

# docker image to use
DOCKER_IMAGE=vryy/kratos_bcn2-cutiga-dev:1.3

# Python command in the container
PY_COMMAND="python3.11"

if [ "$#" -eq 0 ] ; then
    echo "Run all tests"
    export OUTPUT=./ztest_logs/$(whoami)_$(hostname)_docker_"$(date +%Y_%m_%d)".log
    docker run -t -i -v ${HOST_DIR}:/work -w=/work \
        -e OMP_NUM_THREADS=${NUM_THREADS} \
        -v kratos_libs:/home/kuser/kratos_libs \
        -v ${HOME}/workspace/docker/docker-cutiga-dev/kratos_bcn2:/home/kuser/kratos_bcn2 \
        ${DOCKER_IMAGE} /bin/sh -c ". /opt/intel/oneapi/setvars.sh; cd ${SIMULATION_DIR}; python3 run_tests.py $PY_COMMAND | tee $OUTPUT"
elif [ "$#" -eq 1 ] ; then
    echo "Run tests with tag" $1
    export OUTPUT=./ztest_logs/$(whoami)_$(hostname)_docker_"$(date +%Y_%m_%d)"-$1.log
    docker run -t -i -v ${HOST_DIR}:/work -w=/work \
        -e OMP_NUM_THREADS=${NUM_THREADS} \
        -v kratos_libs:/home/kuser/kratos_libs \
        -v ${HOME}/workspace/docker/docker-cutiga-dev/kratos_bcn2:/home/kuser/kratos_bcn2 \
        ${DOCKER_IMAGE} /bin/sh -c ". /opt/intel/oneapi/setvars.sh; cd ${SIMULATION_DIR}; python3 run_tests.py $PY_COMMAND $1 | tee $OUTPUT"
fi
