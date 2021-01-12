#!/bin/sh
BASENAME=`basename "$0"`
USAGE="USAGE: $BASENAME blant.exe M TOP train.el test.el [list of k values, default k=4]"
die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }

[ $# -ge 5 ] || die "need at least 5 args"
#[ -x "$1" ] || die "first arg must be executable"
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
	hawk '{ printf "%s\t%s\n", MAX($1,$2),MIN($1,$2) }' |
	fgrep -c -f - "$test"
) | awk 'NR==1{K=$0}
	    NR>1{printf "\t\t\t**********************\tM='$M' K=%s: Recovered %d edges out of '$TOP', precision %g%% ********************\n", K, $1,100*$1/'$TOP'}'
