#!/bin/bash

export ENV_ROOT="/nvme/marpin/MGRB_phase2_hs37d5x"
export TMPDIR=/nvme/tmp

source "${ENV_ROOT}/software/py2env/bin/activate"

export SPARK_HOME="${ENV_ROOT}/software/spark-2.1-bin-hadoop2.7"
export HAIL_HOME="${ENV_ROOT}/software/hail-0.1-spark-2.1"
export PYTHONPATH="$PYTHONPATH:$HAIL_HOME/python:$SPARK_HOME/python:`echo $SPARK_HOME/python/lib/py4j*-src.zip`"
export SPARK_EXECUTOR_MEMORY=17G
export PYSPARK_SUBMIT_ARGS="--driver-memory 500g --driver-class-path $HAIL_HOME/build/libs/hail-all-spark.jar --jars \"${ENV_ROOT}/software/slf4j-1.7.25/slf4j-simple-1.7.25.jar\" --conf spark.network.timeout=600s --conf spark.executor.cores=1 pyspark-shell"
export SPARK_LOCAL_DIRS="${TMPDIR}/spark"

python "$@"
