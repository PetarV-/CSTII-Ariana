#!/bin/bash

echo "Testing NSGA-II."

mkdir test_results &> /dev/null
tests=("SCH" "FON" "POL" "KUR" "ZDT1" "ZDT2" "ZDT3" "ZDT4" "ZDT6")
testobj_dir='test_objs'
testin_dir='test_inputs'
testres_dir='test_results'
bin_dir='../../bin'

for test_name in "${tests[@]}"
do
	echo "Testing on $test_name..."
	make clean &> /dev/null
	make OBJFILE="$testobj_dir/objectives_$test_name.cpp" &> /dev/null
$bin_dir/nsga2_tester "$testin_dir/sample_input_$test_name.in" "$testres_dir/result_$test_name.out"
	echo "Testing completed. Output written to test_results/result_$test_name.out"
done

make clean &> /dev/null