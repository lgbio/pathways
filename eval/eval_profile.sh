# Profile for the new protein evaluation framework
export EVAL_HOME=/home/lg/cloud/eval

# Used for local-rmsd metric
EVAL_MMTSB=$EVAL_HOME/tools/programs/mmtsb_toolset
# Used for rigidity Analysis (Rigidity Analysis-degrees of freedom)
export EVAL_FIRST=$EVAL_HOME/tools/programs/first
# Used for  Hydrogen bonds naccess program
export EVAL_NACCESS=$EVAL_HOME/tools/program/naccess
# Used for structural similarity
export EVAL_ALPHABET=$EVAL_HOME/tools/metrics/structural_similarity/alphabet.sa

# Set binaries to PATH
EVAL_BIN=$EVAL_HOME/tools/programs/bin:$EVAL_HOME/tools/metrics/bin
EVAL_BIN=$EVAL_BIN:$EVAL_MMTSB:$EVAL_FIRST:$EVAL_NACCESS

export PATH=$PATH:$EVAL_BIN
