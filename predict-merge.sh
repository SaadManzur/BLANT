#!/bin/sh
USAGE="USAGE: $0 {list of blant -mp output files} [DEPRECATED, USE './blant -mq -kK' instead, where K is the largest input k]
PURPOSE: given a bunch of 'blant-mp' output files, merge them so that the node pair u:v is sorted so that u<v,
    and for each canonical-node-pair 'cnp' associated with u:v, the counts C[uv][cnp] are all merged together.
    Both the input and the output files have the following format:
u:v e	k:g:i:j c1	k:g:i:j c2	k:g:i:j c3	[etc]
    Note that we assume TAB-separated columns of (k:o:p c) pairs, and within the column the k:o:p and its count
    are separated by a SPACE. (So you need to distinguish between tabs and spaces on each line.)"

hawk -T '{
	WARN(NF>1,"expecting more than one column");
	# separate u:v from edge value:
	WARN(2==split($1,b," "),"expecting space-separated (u:v edge) value at first column");
	# separate u from v
	WARN(2==split(b[1],a,":"), "expecting colon-separated u:v pair as first element of line");
	if(a[1]<a[2]) uv=b[1]; else uv=a[2]":"a[1]
	if(uv in e) ASSERT(e[uv]==b[2],"node-pair "uv" has different edge truth value than previously seen");
	else e[uv]=b[2];
	for(i=2;i<=NF;i++) {
	    # separate canonical node pair from count
	    WARN(2==split($i,b," "),"expecting space-separated (canonical node-pair count) value at column "i);
	    # The following split is just for syntax checking.
	    WARN(4==split(b[1],a,":"), "expecting colon-separated k:g:i:j quadruplet at column "i);
	    cnp[uv][b[1]]+=b[2];
	}
    }
    END{
	for(uv in e) {
	    printf "%s %d",uv,e[uv];
	    for(cnpc in cnp[uv]) printf "\t%s %d",cnpc,cnp[uv][cnpc]
	    print ""
	}
    }' "$@" || echo "$USAGE" >&2
