#!/bin/bash -ue
PCG -str generate.str -Bondlength  0.2  -LLIB /data1/jackson/TS2CG_pipeline/top/Martini3+NLs.LIB -function  analytical_shape
rm pcg.log

cat header.txt > system.top
cat output.top >> system.top

sed -i 's/PATH//data1/jackson/TS2CG_pipeline/g' system.top
