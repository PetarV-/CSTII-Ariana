#!/bin/bash

echo "Generating synthetic data."

params_dir='parameters'

make clean &> /dev/null
make syn_gen &> /dev/null
../../bin/syn_gen $1 "$params_dir/parameters_$2" "training.out" &> /dev/null
../../bin/syn_gen $1 "$params_dir/parameters_$2" "test.out" &> /dev/null

echo "Synthetic data generated and written to training.out and test.out."

make clean &> /dev/null