#!/bin/bash

echo "Generating synthetic data."

params_dir='parameters'

make clean &> /dev/null
make syn_gen
./syn_gen "$params_dir/$1" &> /dev/null
mv synthetic.out training.out
./syn_gen "$params_dir/$1" &> /dev/null
mv synthetic.out test.out

echo "Synthetic data generated and written to training.out and test.out."

make clean &> /dev/null