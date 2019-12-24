#!/bin/bash
for i in 0 1 2 3 
do 
	for e in 1e-7 1e-14 
	do 
		./a.out 4000 $e 
		echo $i 
	done 
done
