#!/bin/bash -x

gmx pdb2gmx -v -f ../1_attach_Dye/setup.pdb -noter -ff amber03ws_dyes -water tip4p
perl top2psf.pl -p topol.top -o topol.psf

