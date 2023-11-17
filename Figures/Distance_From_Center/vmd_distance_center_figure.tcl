display resetview
color Display Background white
axes location Off
display shadows on
display ambientocclusion on
display aodirect 0.400000
display projection Orthographic
scale by 1.500000
light 0 off
display resize 1000 1000

source /data1/jackson/MD/forcefields/martini3FF/cg_bonds-v2020.tcl
cg_bonds -top modified.top

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
mol modselect 0 0 z < 280
mol modstyle 0 0 Licorice 1.500000 12.000000 12.000000
mol modcolor 0 0 ColorID 2
mol modmaterial 0 0 Smooth

mol addrep 0
mol modselect 1 0 resid 671
mol modstyle 1 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 1 0 ColorID 10
mol modmaterial 1 0 Smooth

mol addrep 0
mol modselect 2 0 resid 4743
mol modstyle 2 0 Licorice 2.000000 12.000000 12.000000
mol modcolor 2 0 ColorID 3
mol modmaterial 2 0 Smooth

render TachyonLOptiXInternal top.bmp -res 4000 4000 -format bmp

quit