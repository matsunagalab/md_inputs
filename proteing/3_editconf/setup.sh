#!/bin/bash -x

#gmx editconf -f ../2_pdb2gmx/conf.gro -c -d 2.00 -bt cubic
gmx editconf -f ../2_pdb2gmx/conf.gro -c -box 10.50 -bt cubic

