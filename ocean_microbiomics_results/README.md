# Ocean Microbiomics

## The dataset

The dataset we are working with is the 34,000+ genomes (refs, sags, and mags) from the ocean microbiomics database v1.

## The scripts

There are 2 scripts in this repo. The processingresults.R takes an input file - the hmm search results - and determines the steps present for each pathway in each genome. The processingresults_new_filtered.R takes in the new results that use new hmm score cutoffs (less stringent for sulfate assimilation) and redoes the calculations and figures for sulfate assimilation.

## Sulfate.xlsx

This file contains the information for each pathway necessary to determine completion percentages and such.

## The directories

The two directories contain the files produced in the two scirpts above. There are csv files (raw and processed) and a myriad of pdfs that help to visualize the results.
