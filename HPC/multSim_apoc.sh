#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=18G
#$ -l h_rt=240:0:0
#$ -wd ~/HSCDynamics
#$ -j y
#$ -N somaticMoranMultSimulations
#$ -o ~/HSCDynamics
#$ -t 1-99

module load julia

Nf=1000

SRCDIR=$HOME/HSCDynamics/
DATADIR=$HOME/HSCDynamics/data/Nf$Nf
mkdir -p $DATADIR

rsync -rltv $SRCDIR/src $TMPDIR/
rsync -rltv $SRCDIR/HPC $TMPDIR/

cd $TMPDIR

echo "starting julia script..."

# julia HPC/addPackages.jl
julia HPC/singleSimScript.jl ${SGE_TASK_ID} $Nf

# rsync -rltv $TMPDIR/ $DATADIR/
mv singlePatientFullSim_*.jld2 $DATADIR/

echo "success"

