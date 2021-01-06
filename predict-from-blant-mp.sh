#!/bin/sh
USAGE="$0 [-include-known] blant-mp-file
PURPOSE: given a blant-mp output file, learn which motifs have predictive value, and then use the precision curves to
create a list of predictions sorted best-to-worst. By default, we only output predictions on the set of node pairs
that had *no* edge in the blast-mp file; these are genuine predictions. If the '-include-known' option is given, then
the 'prediction' is included even if the edge was already in the input data. This facilitates measuring precision on
the training data."

die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }

TEST="$1"
INCLUDE_KNOWN=0
while [ $# -gt 1 ]; do
    case "$1" in
    -include-known) INCLUDE_KNOWN=1; shift;;
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

[ $# -eq 1 ] || die "expecting a blant-mp output file"

# input lines look like:
#ENSG00000197362:ENSG00000204178 0	4:9:0:4 4	4:5:1:2 46	4:6:0:3 6	4:7:1:2 2	[ etc ... ]
# where the colon word is k:g:i:j (g=graphlet Ordinal, i,j is a node (NOT ORBIT) pair in g.), followed by a count.
# We call (i,j) a "canonical node pair", or cnp for short.

hawk 'BEGIN{min_samples=1000; min_rho=0.6; min_t=100; min_p=0.9} # THESE MAY NEED TO BE ADJUSTED
    ARGIND==1{
	uv=$1 # node pair
	E[uv]=e=$2 # edge Boolean
	ASSERT($2==0 || $2==1, "expecting second column to be Boolean");
	for(i=3;i<NF;i+=2){ #col 3 onwards are (cnp,count) pairs
	    cnp=$i; c=$(i+1)
	    PearsonAddSample(cnp,c,e);
	    ++hist[cnp][c][e]; # histogram: number of times orbit-pair $i had exactly count c for edge truth "e"
	    if(c>max[cnp])max[cnp]=c
	}
    }
    ENDFILE{if(ARGIND==1){
	# Now produce empirical precision curve as a function of orbit-pair count (c), for each orbit-pair
	printf("Predictive stats for %s\n", FILENAME) > "/dev/stderr"
	for(cnp in hist) {
	    for(c=max[cnp];c>=0; --c) { # starting at the highest orbit-pair counts
		numer[cnp][c]+=hist[cnp][c][1]; # those that actually had edge
		denom[cnp][c]+=hist[cnp][c][1]+hist[cnp][c][0]; # both edge and non-edge
		# precision of this orbit-pair as func of c (add 1 to denom simply to avoid div by zero)
		prec[cnp][c] = numer[cnp][c]/(denom[cnp][c]+1);
	    }
	    # Now, a heuristic "fix" to statistical noise for the few, highest-scoring counts: fix the top predictions
	    # so that the precision is artificially fixed to be monotonically increasing with orbit pair count
	    maxP[cnp]=0;
	    for(c=0;c<=max[cnp];c++) {
		# make precision non-decreasing with increasing orbit count
		prec[cnp][c]=MAX(maxP[cnp],prec[cnp][c]);
		maxP[cnp]=MAX(maxP[cnp],prec[cnp][c]);
	    }
	    if(_Pearson_N[cnp]>min_samples){
		PearsonCompute(cnp);
		if(maxP[cnp] > min_p && _Pearson_rho[cnp]>min_rho && _Pearson_t[cnp] > min_t)
		    print PearsonPrint(cnp),cnp,maxP[cnp] > "/dev/stderr"
	    }
	}
    }}
    ARGIND==2{ # actually the same file, just that we go through it now creating predictions
	uv=$1 # node pair
	p1=0; e=$2;
	c=0; bestCol="none";
	for(i=3;i<NF;i+=2){
	    cnp=$i; c=$(i+1);
	    if(c>max[cnp])c=max[cnp]; # clip the orbit count to the highest seen during training
	    # filter on "reasonable" orbit pairs
	    if(_Pearson_N[cnp]>min_samples && _Pearson_rho[cnp] > min_rho && _Pearson_t[cnp] > min_t && maxP[cnp]>=min_p)
		if (p1<prec[cnp][c]){p1=prec[cnp][c];bestCol=c"="cnp} #... and take the best resulting prediction
	}
	if(p1>min_p && E[uv]<='$INCLUDE_KNOWN') printf "%s\t%g\tbestCol %s\t[%s]\n",uv,p1,bestCol,$0
    }' "$1" "$1" | # yes, twice
    sort -S4G -k 2gr -k 4nr
    #awk '{sum+=$6; printf "%d prec %g %s\n", NR, sum/NR,$0}'
