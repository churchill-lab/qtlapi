#!/usr/bin/env bash
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=6g
#SBATCH --time=02:00:00

module load R/3.6.0

#cd /fastscratch/mvincent/

echo "FINDPEAKSR : $FINDPEAKSR"
echo "DATAFILE : $DATAFILE"
echo "QTLAPISOURCE : $QTLAPISOURCE"
echo "OUTPUTFILE : $OUTPUTFILE"
echo "DATASET : $DATASET"
echo "START : $START"
echo "STEP : $STEP"

Rscript $FINDPEAKSR $DATAFILE $QTLAPISOURCE $OUTPUTFILE $DATASET $START $STEP


