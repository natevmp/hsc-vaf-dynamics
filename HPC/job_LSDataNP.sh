#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -l h_rt=48:0:0
#$ -wd ~/HSCDynamics
#$ -j y
#$ -N LSDataNP
#$ -o ~/HSCDynamics
#$ -m beas
#$ -t 1-3

module load julia

SRCDIR=$HOME/HSCDynamics/
DATADIR=$HOME/HSCDynamics/data/
mkdir -p $DATADIR

rsync -rltv $SRCDIR/src $TMPDIR/
rsync -rltv $SRCDIR/HPC $TMPDIR/

cd $TMPDIR

echo "starting julia script..."

# julia HPC/addPackages.jl
julia HPC/LSData_N-p_Space.jl ${SGE_TASK_ID}

# rsync -rltv $TMPDIR/ $DATADIR/
mv LSData_NPSpaceInference_tM*.jld2 $DATADIR/

echo "success"

