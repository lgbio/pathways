#!/bin/bash

# Calculates the Potential energy for a protein conformation
# It uses gromacs g_energy and it needs to construct a full trajectory files

PDB_REFERENCE_FILE=$1
PDB_PROTEIN_FILE=$2

REFERENCE=`echo $PDB_REFERENCE_FILE | cut -d "." -f 1` 
PROTEIN=`echo $PDB_PROTEIN_FILE | cut -d "." -f 1` 
echo ">>>" $REFERENCE $PROTEIN  >& $PROTEIN.log

EXTRA_OPTIONS='-ignh'

pdb2gmx -f $PROTEIN.pdb -p $PROTEIN.top -o $PROTEIN.gro -i $PROTEIN.itp -ff G43a1 $EXTRA_OPTIONS >> $PROTEIN.log 2>&1

touch $PROTEIN.mdp 

grompp  -f $PROTEIN.mdp -po /dev/null -o $PROTEIN.tpr -c $PROTEIN.gro -p $PROTEIN.top >> $PROTEIN.log 2>&1

mdrun  -s $PROTEIN.tpr -o $PROTEIN.trr -c /dev/null -e $PROTEIN.edr -cpo $PROTEIN.cpt -g /dev/null>> $PROTEIN.log 2>&1

g_mindist -f $PROTEIN.trr -s $REFERENCE.gro -ox $PROTEIN.trr -d 0.7 -o $PROTEIN.xvg
#g_mindist -f $PROTEIN-min.trr -s $REFERENCE.gro -o $PROTEIN.xvg -d 0.3 -max

echo `tail -1 $PROTEIN.xvg`|split-lg.py "" -1

