#!/bin/bash -x

gmx pdb2gmx -v -f ../1_attach_Dye/setup.pdb -noter -ff amber99sb-ildn_dyes -water tip3p -o conf.pdb
perl top2psf.pl -p topol.top -o topol.psf

