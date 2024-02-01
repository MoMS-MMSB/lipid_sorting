#!/bin/bash -ue
PCG -str generate.str -Bondlength  0.2  -LLIB /data2/jackson/lipid_sorting/lipid_sorting_in_tubules/TS2CG-Setup-Pipeline/top/Martini3+NLs.LIB -function  analytical_shape
rm pcg.log

cat header.txt > system.top
cat output.top >> system.top
