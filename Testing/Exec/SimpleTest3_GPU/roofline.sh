#!/bin/bash

# This script grabs DRAM bytes read+write and all floating point ops so that
# you can calculate arithmetic intensity for a roofline plot.

srun -n 1 -c 2 --cpu-bind=cores \
  nv-nsight-cu-cli \
  --metrics dram__bytes_read.sum,dram__bytes_write.sum,smsp__sass_thread_inst_executed_op_dadd_pred_on.sum,smsp__sass_thread_inst_executed_op_dmul_pred_on.sum,smsp__sass_thread_inst_executed_op_dfma_pred_on.sum \
  Pele3d.gnu.TPROF.CUDA.ex \
  inputs.2d
