#!/bin/bash
#SBATCH --array=1-121
# #SBATCH --array=1-1
#SBATCH -o logs/scan_%A_%a.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH --overcommit
#SBATCH --ntasks-per-core=2
#SBATCH -t 119:59:00
#SBATCH --mem=4000

OFFSET=1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)
rundir='full__1a__model_1__log10c0_-2to2__p1_0.5to1__E_1__ctrl0_1__log10(delta_to_c0)_-2'

outdir="../Data/Raw/${rundir}"
line=$(sed -n "$LINE_NUM"p $outdir/to_run.csv)

# log10c0,p1,E,log10delta

log10c0=$(echo "$line" | cut -d "," -f 1)
p1=$(echo "$line" | cut -d "," -f 2)
E=$(echo "$line" | cut -d "," -f 3)
log10delta=$(echo "$line" | cut -d "," -f 4)
outfile="$outdir/out_${LINE_NUM}.mat"

echo "Offset $OFFSET ; Line $LINE_NUM ; log10c0 $log10c0 ; p1 $p1 ; E $E ; log10delta $log10delta ; outfile $outfile"

matlab -r "load('$outdir/params.mat'); \
           params.log10c0 = $log10c0; \
           params.P = [$p1; 1-$p1]; \
           params.E = $E'; \
           params.log10delta = $log10delta'; \
           disp('A serial-dilution simulation is starting...'); \
           output = sim__serial__full_dynamics(params,'$outfile'); \
           disp('Finished!'); \
           quit();"
