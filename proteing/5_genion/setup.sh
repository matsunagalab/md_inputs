#!/bin/bash -x

cp ../4_solvate/topol.top .
gmx grompp -f ions.mdp -c ../4_solvate/out.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -np 8

