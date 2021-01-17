#!/bin/sh
USAGE="$0 k
PURPOSE: extract which orbits are predictive for HI-union from previous runs in ~BLANT/HI-union"
die(){ (echo "USAGE: $USAGE"; echo "`basename $0`: FATAL ERROR: $@")>&2; exit 1; }

k=$1; shift
[ 4 -le "$k" -a "$k" -le 8 ] || die "first argument must be k, between 4 and 8"

awk 'BEGIN{k='$k'
	min_n=10000; min_rho=0.1; min_s=60; min_p=.7
	#if(k==7)     min_rho=0.3  # to reduce RAM usage to tolerable levels :-(
	}
    /:/{
	n=1*$1; rho=1*$2; pv=1*$3; s=1*$4; p=1*$6;
	if(n>min_n && rho>min_rho && s>min_s && p>min_p)print $5
    }' /home/wayne/src/bionets/BLANT/HI-union/*k$k*.predictors.loose | sort -u
