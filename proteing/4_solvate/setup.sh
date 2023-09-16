#!/bin/bash -x

cp ../2_pdb2gmx/topol.top .
gmx solvate -cp ../3_editconf/out.gro -cs amber03ws_dyes.ff/tip4p2005.gro -p topol.top

