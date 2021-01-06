#!/bin/sh
USAGE="USAGE: $0 blant.exe M TOP train.el test.el test.el [list of k values]"

die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $@}" </dev/null; }

NL='
'
TAB='	'

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

[ -x "$1" ] || die "first arg must be executable"
BLANT="$1"
M=$2
TOP=$3
train=$4
test=$5
shift 5
( echo "$@" # just for awk to print out later
    ./predict-edges-from-network.sh -B "$BLANT" -M $M "$train" "$@" |
	sed 's/:/    /' |
	head -$TOP |
	awk '{ printf "%s\t%s\n%s\t%s\n", $1,$2,$2,$1 }' |
	fgrep -c -f - "$test"
) | awk 'NR==1{K=$0}
	    NR>1{printf "\t\t\t**********************\tM='$M' K=%s: Recovered %d edges out of '$TOP', precision %g%% ********************\n", K, $1,100*$1/'$TOP'}'
