mkdir 1.initial
cd 1.initial
nextflow run ../../TS2CG-Setup-Pipeline/main.nf 

cd ../
mkdir 2.eq_nopore
gmx grompp -f eq-nopore.mdp -c ../1.initial/results/eq.gro -p ../1.initial/results/topol.top -n ../1.initial/results/index.ndx -o eq_nopore.tpr

# gmx mdrun -v -deffnm eq-nopore