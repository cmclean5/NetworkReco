#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/NetworkReco/v1r1
SCRIPTDIR=/exports/home/cmclean5/STUDIES/NetworkReco/v1r1/EDDIE/SCRIPTS

N=1
RESTARTS=10

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

cp -r $SCRIPTDIR/EM.R .
cp -r $SCRIPTDIR/Generate_Rndm_Measurements.R .
cp -r $SCRIPTDIR/test_monomodal.R .

# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
module load R 
#module load igmm/apps/R/3.2.2 #this should work for igraph & poweRlaw

time Rscript test_monomodal.R $JOB_ID $N $RESTARTS

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.csv $OUTDIR

echo "$0 done!"

exit 0
