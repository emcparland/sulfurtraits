#!/bin/sh

n=0

awk 'NR>1 {print $1}' customhmm_scorecutoffs.txt | while read headername

do
	n=n+1
	score=$(awk -v var="$headername" '$1 == var' customhmm_scorecutoffs.txt | awk '{print $4}')
	cat $headername | awk -v var=$score 'NR>3 && $6 > var' | sed '/^#/d' >> allcustom_filtered.txt
done
echo $n
