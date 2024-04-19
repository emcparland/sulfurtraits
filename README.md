# sulfurtraits

##definition of scripts
### kofamsearch
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
### kofamsearch custom
kofamsearch_custom.sh: same format as above, except this is for some custom HMM profiles that we created for n = 22 genes that are not in KO but we are interestedin.

filter_customhmm.sh: will filter the output with our custom score values and create one concatanated file

Checking custom evaluves and scores nb: plots in R to look at the chosen score value cutoffs (using third quartile)

customhmm_scorecutoffs.txt: contains names of files, gene of interest, the original cutoff value, and the more stringent cutoff value (which == the 3rd quartile)

###data results
All stored in our google drive for now under data_for_manuscript

all_filtered_results.txt: the filtered results for the 240 KOs from Kegg as described above.

allcustom_filtered.txt: the filtered results for the custom HMMs as described above.
