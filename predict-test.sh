#!/bin/sh
TOP=500
M=8
while true;
do
    for K in 4 '4 5' '4 5 6' '4 5 6 7' '4 5 6 7';
    do
	time ./predict-edges-from-network.sh -B ./blant.test -M $M ~/extra1/preserve/Kovacs/HI-union/split-by-edges/F0/HI-union-train0.el $K |
	    sed 's/:/    /' |
	    head -$TOP |
	    hawk '{ printf "%s\t%s\n", MAX($1,$2),MIN($1,$2)}' |
	    fgrep -c -f - ~/extra1/preserve/Kovacs/HI-union/split-by-edges/F0/HI-union-test0.el |
	    awk '{printf "\t\t\t**********************\tM='$M' K='"$K"': Recovered %d edges out of '$TOP', precision %g%% ********************\n",
		$1,100*$1/'$TOP'}'
    done
    M=`expr $M '*' 2`
done
