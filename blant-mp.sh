#!/bin/sh
USAGE="USAGE: $0 [-async] -C {2|3|4|5} blant-exec {then blant args, including -mp}
[2345] is COLS, and mandatory: which columns to remember?"
die(){(echo "$USAGE"; echo "$0: FATAL ERROR: $@") >&2; exit 1; }

ASYNC=false
while [ $# -gt 0 ]; do
    case "$1" in
    -a*) ASYNC=true; shift;;
    -C) COLS=$2; shift 2;;
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

echo "$COLS" | grep -q '^[2345]$' || die "expecting COLS to be one of [2345]"

k=`echo "$@" | sed -e 's/.*-k//' | awk '{printf "%d\n",$1}'`
echo "$k" | grep -q '^[45678]$' || die "expecting -k argument with 4-8"

BLANT=$1; shift

$BLANT "$@" | hawk 'BEGIN{C='$COLS'}
    {
	#blant raw output looks like: (leading P just because...)
	#P  u:v    e  o:p   q:r   x:y
	#P 680:265 0 37:37 38:39 517:476
	#P 680:265 0 37:37 38:39 517:25
	#P 265:517 1 37:38 37:39 680:476
	sub("^","'$k':",$4);
	uv=$2; e=$3; op=$4; qr=$5; xy=$6
	G[uv]=e;
	# PG=PredictGraph, but what to store and count? maybe q:r is irrelevant, just unique edges?
	if(C==5)        PG[uv][op][qr][xy]=1 # worst in terms of memory usage: distinguish the whole (u,v,o,p,q,r,x,y) octuplet
	else if(C==4)   PG[uv][op][qr]=1 # this one ignores the (x,y) edge and only counts (q,r) pairs
	else if(C==3)   PG[uv][op][xy]=1 # this one ignores the (q,r) orbit and only counts nearby edges.
	else if(C==2) ++PG[uv][op]  # Simplest: literally just count (u:v,o:p) orbit-pair frequencies
	else ASSERT(false,"wrong value of C???");
    }
    END {
	for(uv in PG) {
	    printf "%s %d",uv,G[uv]
	    for(op in PG[uv]) {
		if(C==5) {
		    value = 0 # TP=TotalPredict of edge-orbit-pair (o,p) on node pair (u,v) in G
		    for(qr in PG[uv][op]) value += length(PG[uv][op][qr]) # all the 1s of the xys
		}
		else if(C==4||C==3) value = length(PG[uv][op])
		else if(C==2) value = PG[uv][op]
		printf "\t%s %d",op,value;
	    }
	    print ""
	}
    }' &
if $ASYNC; then
    pid=`jobs -l | head -1 | awk '{print $2}'`

    while free -g | awk '/Mem:/&&$NF<2{print;++flag}/Swap:/&&$NF<2{print;++flag}flag==2{exit(flag)}'; do
	sleep 10; free -g|head -1;
    done; kill $pid # nuke the $BLANT generating the stuff for hawk
fi
wait
