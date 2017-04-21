#!/bin/bash

# Calculates the RMSD value for two protein structures
# It uses gromacs g_rms and it needs to construct gromacs files

# NOTES:
#	- Fails when some PDBs are missing atoms

PDB_REFERENCE_FILE=$1
PDB_INPUT_FILE=$2

PROTEIN=`echo $PDB_INPUT_FILE | cut -d "." -f 1` 
REFERENCE=`echo $PDB_REFERENCE_FILE | cut -d "." -f 1` 

echo ">>>" g_rms $REFERENCE $PROTEIN >& $PROTEIN.log

EXTRA_OPTIONS='-ignh'

case $PHD_FORCE_FIELD in
	"AMBER96") FORCE_FIELD=amber96;;
	"OPLSAA")  FORCE_FIELD=oplsaa;;
	"GROMOS")  FORCE_FIELD=G43a1;;
	*) echo "ERROR in rmsd: Invalid Force Field"; exit;;
esac

echo "pdb2gmx -f $REFERENCE.pdb -p $REFERENCE.top -o $REFERENCE.gro -i $REFERENCE.itp -ff $FORCE_FIELD $EXTRA_OPTIONS >> $PROTEIN.log 2>&1"
pdb2gmx -f $REFERENCE.pdb -p $REFERENCE.top -o $REFERENCE.gro -i $REFERENCE.itp -ff $FORCE_FIELD $EXTRA_OPTIONS >> $PROTEIN.log 2>&1
echo "pdb2gmx -f $PROTEIN.pdb -p $PROTEIN.top -o $PROTEIN.gro -i $PROTEIN.itp -ff $FORCE_FIELD $EXTRA_OPTIONS >> $PROTEIN.log 2>&1"
pdb2gmx -f $PROTEIN.pdb -p $PROTEIN.top -o $PROTEIN.gro -i $PROTEIN.itp -ff $FORCE_FIELD $EXTRA_OPTIONS >> $PROTEIN.log 2>&1

echo C-alpha C-alpha | g_rms -s $REFERENCE.gro -f $PROTEIN.gro -o $PROTEIN.xvg  >> $PROTEIN.log 2>&1

lastline=`tail -1 $PROTEIN.xvg`

strvalue=`split-lg.py "" -1 "$lastline"`

amnstrongValue=`calc.py $strvalue 10.0`

echo $amnstrongValue

rm $PROTEIN.gro $PROTEIN.itp $PROTEIN.log $PROTEIN.top $PROTEIN.xvg
rm $REFERENCE.gro $REFERENCE.itp $REFERENCE.top

