for dir in */; do 
   cd $dir 
   vmd *gro -dispdev none -startup ../vmd_tubule_figures.tcl 
   cd ../
done