#!/bin/bash

echo "Testing NSGA-II."

mkdir test_results &> /dev/null
tests=("SCH" "FON" "POL" "KUR" "ZDT1" "ZDT2" "ZDT3" "ZDT4" "ZDT6")
src_dir='../../src/nsga2'
testobj_dir='../../test/nsga2/test_objs'
testin_dir='../../test/nsga2/test_inputs'
testres_dir='../../test/nsga2/test_results'

cd $src_dir

for test_name in "${tests[@]}"
do
	echo "Testing on $test_name..."
	make clean &> /dev/null
	make OBJFILE=$testobj_dir/objectives_$test_name.cpp &> /dev/null
	./nsga2 "$testin_dir/sample_input_$test_name" &> /dev/null
	mv results.out $testres_dir/result_$test_name.out &> /dev/null
	echo "Testing completed. Output written to test_results/result_$test_name.out"
done

make clean &> /dev/null