#!/bin/bash

# Calculates the Potential energy for a tmpFile conformation
# It uses gromacs g_energy and it needs to construct a full trajectory files

INPUT_FILE=$1
PROTEIN=`echo $INPUT_FILE | cut -d "." -f 1` 
tmpFile=$PROTEIN-energy-tmp
echo ">>> energy-lg.sh " $tmpFile >& $tmpFile.log

EXTRA_OPTIONS='-ignh'

pdb2gmx -f $PROTEIN.pdb -p $tmpFile.top -o $tmpFile.gro -i $tmpFile.itp -ff G43a1 $EXTRA_OPTIONS >> $tmpFile.log 2>&1

touch $tmpFile.mdp 

grompp  -f $tmpFile.mdp -po $tmpFile.mdp -o $tmpFile.tpr -c $tmpFile.gro -p $tmpFile.top >> $tmpFile.log 2>&1

mdrun  -s $tmpFile.tpr -o $tmpFile.trr -c $tmpFile.gro -e $tmpFile.edr -cpo $tmpFile.cpt -g $tmpFile.log>> $tmpFile.log 2>&1

echo Potential | g_energy -f $tmpFile.edr -o $tmpFile.xvg >> $tmpFile.log 2>&1

if [ -f "$tmpFile.xvg" ]
then
	echo `tail -1 $tmpFile.xvg`|split-lg.py "" -1
else
	echo 9999999999999999.9999
fi

echo ">>>" $tmpFile
rm *$tmpFile*
