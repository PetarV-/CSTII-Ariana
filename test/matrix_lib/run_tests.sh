#!/bin/bash

echo "Testing Matrix Exponential."

mkdir test_results &> /dev/null
tests=("Moler1" "MolerVanLoan1" "MolerVanLoan2" "MolerVanLoan3" "MolerVanLoan4" "Ward1" "Ward2" "Ward3" "Ward4")
testin_dir='exp_tests'
testres_dir='test_results'

make

for test_name in "${tests[@]}"
do
	echo "Testing on $test_name..."
	./exp_tester "$testin_dir/input_$test_name" &> /dev/null
	mv exp_result.out $testres_dir/result_$test_name.out &> /dev/null
	echo "Testing completed. Output written to test_results/result_$test_name.out"
done

make clean &> /dev/null