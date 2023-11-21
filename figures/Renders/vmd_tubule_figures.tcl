display resetview
color Display Background white
axes location Off
display shadows on
display ambientocclusion on
display aodirect 0.400000
display projection Orthographic
scale by 1.500000
light 0 off
# display dof on
display resize 1000 1000

source /data1/jackson/MD/forcefields/martini3FF/cg_bonds-v2020.tcl
cg_bonds -top topol.top

material add copy RTChrome
material rename Material23 Smooth
material change ambient Smooth 0.020000
material change specular Smooth 0.200000
material change mirror Smooth 0.200000
material change shininess Smooth 0.200000
material change outline Smooth 0.000000
material change outlinewidth Smooth 0.000000

mol delrep 0 0

mol addrep 0
mol modselect 0 0 resname DOPC
mol modstyle 0 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 0 0 ColorID 10
mol modmaterial 0 0 Smooth

mol addrep 0
mol modselect 1 0 resname POPC
mol modstyle 1 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 1 0 ColorID 15
mol modmaterial 1 0 Smooth

mol addrep 0
mol modselect 2 0 resname POPE
mol modstyle 2 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 2 0 ColorID 3
mol modmaterial 2 0 Smooth

mol addrep 0
mol modselect 3 0 resname DOPE
mol modstyle 3 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 3 0 ColorID 17
mol modmaterial 3 0 Smooth

mol addrep 0
mol modselect 4 0 resname DOPS
mol modstyle 4 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 4 0 ColorID 28
mol modmaterial 4 0 Smooth

mol addrep 0
mol modselect 5 0 resname DOPA
mol modstyle 5 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 5 0 ColorID 25
mol modmaterial 5 0 Smooth

# display dof_focaldist 0.50000
render TachyonLOptiXInternal top.bmp -res 2000 2000 -format bmp

rotate x by 90.000000
rotate y by 90.000000
render TachyonLOptiXInternal rotated.bmp -res 2000 2000 -format bmp

rotate x by 30.000000
render TachyonLOptiXInternal tilted.bmp -res 2000 2000 -format bmp

quit