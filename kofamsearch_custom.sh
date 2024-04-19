#!/bin/sh
#SBATCH --partition=compute
#SBATCH --job-name=kofam
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=emcparland@whoi.edu
#SBATCH --ntasks-per-node=12
#SBATCH --mem=25gb
#SBATCH --time=24:00:00
#SBATCH --output=kofam_%A-%a.out
#SBATCH --error=kofam_%A-%a.error
#SBATCH --array=1-22%22

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate sulfbio

kofampath=~/kofamscan/kofam_scan-1.3.0
wkdir=~/sulfomics/custom_hmm_mdemers
echo $SLURM_ARRAY_TASK_ID
line1="awk 'NR==$SLURM_ARRAY_TASK_ID {print $"1"}'"
headername=$(cat $wkdir/custom.hal|eval $line1)
echo $headername
line2="$kofampath/./exec_annotation -p $wkdir/profiles/$headername --cpu=12 --tmp-dir=./tmp$SLURM_ARRAY_TASK_ID -f detail-tsv -o $wkdir/$headername.txt genes-redundant.faa"
echo $line2
eval $line2
