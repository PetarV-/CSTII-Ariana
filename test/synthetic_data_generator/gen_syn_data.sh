#!/bin/bash

echo "Generating synthetic data."

params_dir='parameters'

make clean &> /dev/null
make syn_gen &> /dev/null
../../bin/syn_gen "$params_dir/parameters_$1" "training.out" &> /dev/null
../../bin/syn_gen "$params_dir/parameters_$1" "test.out" &> /dev/null

echo "Synthetic data generated and written to training.out and test.out."

make clean &> /dev/null