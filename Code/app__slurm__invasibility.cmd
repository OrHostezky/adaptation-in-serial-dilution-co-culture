#!/bin/bash
# #SBATCH --array=1-81
#SBATCH --array=1-1
#SBATCH -o logs/scan_%A_%a.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH --overcommit
#SBATCH --ntasks-per-core=2
#SBATCH -t 119:59:00
#SBATCH --mem=4000

OFFSET=1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)
rundir='invasibility__1na_1a__model_1__log10c0_-2to2__E_[1_1]__alpha0_[1_0;0.5_0.5]__ctrl0_[1_0]__log10(delta_to_c0)_-3to0__rho_pre_0.99999'

outdir="../Data/Raw/${rundir}"
line=$(sed -n "$LINE_NUM"p $outdir/to_run.csv)

# ctrl0,log10c0,E,alpha0

ctrl0=$(echo "$line" | cut -d "," -f 1)
log10c0=$(echo "$line" | cut -d "," -f 2)
E=$(echo "$line" | cut -d "," -f 3)
alpha0=$(echo "$line" | cut -d "," -f 4)
outfile="$outdir/out_${LINE_NUM}.mat"

echo "Offset $OFFSET ; Line $LINE_NUM ; log10c0 $log10c0 ; E $E ; alpha0 ; ctrl0 $ctrl0 ; outfile $outfile"

matlab -r "load('$outdir/params.mat'); \
           params.log10c0 = $log10c0; \
	   params.E = $E'; \
	   params.alpha0 = $alpha0; \
           params.ctrl0 = $ctrl0; \
	   p1s = 0.5:0.01:1; \
           log10deltas = params.log10c0-flip(0:0.05:3); \
	   rho_pre = 0.99999; \
	   disp('Running invasibility tests...'); \
           isInvaded = sim__invasibility_map(params,'$outfile',p1s,log10deltas,rho_pre); \
           disp('Finished!'); \
           quit();"
