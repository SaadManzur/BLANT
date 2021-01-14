#!/bin/sh
USAGE="$0 k
PURPOSE: extract which orbits are predictive for HI-union from previous runs in ~BLANT/HI-union"
die(){ (echo "USAGE: $USAGE"; echo "`basename $0`: FATAL ERROR: $@")>&2; exit 1; }

k=$1; shift
[ 4 -le "$k" -a "$k" -le 8 ] || die "first argument must be k, between 4 and 8"

awk '/:/{
	n=1*$1;rho=1*$2;pv=1*$3;s=1*$4;p=1*$6;
	if(n>10000 && rho>.1 && pv<1e-9 && s>10 && p>.8)print $5
    }' /home/wayne/src/bionets/BLANT/HI-union/*k$k*.predictors.loose | sort -u
