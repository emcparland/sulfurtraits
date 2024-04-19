# sulfurtraits

##definition of scripts
kofamsearch.sh: will run kofamscan for your KOs of interest using an array, submitting 240 KO's (20 at a time). Important to use the --tmp-dir flag or else the arrays will overwrrite each other and the job fails

when kofamsearch.sh is complete, you can run the following commands to create the concatanated file:

Take advantage of the fact that kofamscan marks lines with a * if they pass the cut-off value.
```
cat hmm_results/*.hmm.txt |wc -l

cat hmm_results/*.hmm.txt | grep "*" > all_filtered_results.txt

cat hmm_results/all_filtered_results.txt |wc -l
```

After filtering, we have eliminated some KOs that were not significant according to kofamscan
```
awk '{print $3}' hmm_results/all_filtered_results.txt |sort |uniq |wc -l
```
(Answer here is 206 KOs)

kofamsearch_custom.sh: same format as above, except this is for some custom HMM profiles that we created
