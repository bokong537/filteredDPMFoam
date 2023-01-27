# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:41:43 2023

@author: bokong
"""

#import os 

px = [-0.005, 0.005]
py = [-0.025, 0.025]
pz = [0, 0.15]

nx = 1
ny = 5
nz = 15

dx = (px[1] - px[0])/(nx+1)
dy = (py[1] - py[0])/(ny+1)
dz = (pz[1] - pz[0])/(nz+1)

file1 = open("kinematicCloudPositions","w")

L = ["/*--------------------------------*- C++ -*----------------------------------*\ \n",
"| =========                 |                                                 |\n",
"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n",
"|  \\    /   O peration     | Version:  v2212                                 |\n",
"|   \\  /    A nd           | Website:  www.openfoam.com                      |\n",
"|    \\/     M anipulation  |                                                 |\n",
"\*---------------------------------------------------------------------------*/\n",
"FoamFile\n",
"{\n",
"    version     2.0;\n",
"    format      ascii;\n",
"    class       vectorField;\n",
"    object      kinematicCloudPositions;\n",
"}\n",
"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
"\n",
"( \n"] 

file1.writelines(L)
 
for i in range(nx) : 
    for j in range(ny):
        for k in range(nz):
            file1.write("(")
            file1.write(str(px[0] + dx*(i+1))+"\t")
            file1.write(str(py[0] + dy*(j+1))+"\t")
            file1.write(str(pz[0] + dz*(k+1))+"\t")
            file1.write(")\n")

file1.write(")\n")
            

file1.write("\n// ************************************************************************* //\n")

file1.close()