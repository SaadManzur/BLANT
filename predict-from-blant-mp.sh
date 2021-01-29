#!/bin/sh
USAGE="$0 [-eval TestEdgeList.el] [-a] [-predictors-only ] [-nonorm] [-include-known] [-sum[log]] [-notL3] blant-mp-file
PURPOSE: given a blant-mp output file, learn which motifs have predictive value, and then use the precision curves to
create a list of predictions sorted best-to-worst. By default, we only output predictions on the set of node pairs
that had *no* edge in the blast-mp file; these are genuine predictions. If the '-include-known' option is given, then
the 'prediction' is included even if the edge was already in the input data. This facilitates measuring precision on
the training data."

die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }

TMPDIR=`mktemp -d /tmp/predict-from-blant-mp.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0
 trap "echo error encountered: TMPDIR is $TMPDIR; exit" 1 2 3 15

WINDOW_SIZE=1
INCLUDE_KNOWN=0
PREDICTORS_ONLY=0
DO_NORM=1 # you normally want to leave this ON because then BINS_PER_WEIGHT_UNIT=100 makes sense; otherwise too many bins
EVALUATE=''
ALL_PREDICTIONS=0
VERBOSE=0
TIGHT_MINIMUMS="min_samples=10000; min_rho=0.15; min_t= 50; min_p=0.8" # When using Pearson
LOOSE_MINIMUMS="min_samples=100  ; min_rho=0.01; min_t=  6; min_p=0.5" # less stringent, just for detection
MINIMUMS="$TIGHT_MINIMUMS" # default is tight

while [ $# -gt 0 ]; do
    case "$1" in
    -include-known) INCLUDE_KNOWN=1; shift;;
    -predictors-only) PREDICTORS_ONLY=1; MINIMUMS="$LOOSE_MINIMUMS"; shift;;
    -v*) VERBOSE=1; shift;;
    -nonorm*) DO_NORM=0; shift;;
    -eval*) EVALUATE="$2"; shift 2;;
    -win*) WINDOW_SIZE="$2"; shift 2;;
    -a*) ALL_PREDICTIONS=1; shift;;
    -) shift; break;; # means we're processing stdin
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

if [ $DO_NORM = 0 ]; then
    THIRD_PASS=''
else
    THIRD_PASS="$TMPDIR/input" # Just the filename; file doesn't exist until it's created below
fi

if [ $# -eq 0 ]; then
    echo processing stdin >&2
    set -
else
    echo processing "$@" >&2
fi

wzcat "$@" > $TMPDIR/input

# input lines look like:
#ENSG00000197362:ENSG00000204178 0	4:9:0:4 4	4:5:1:2 46	4:6:0:3 6	4:7:1:2 2	[ etc ... ]
# where the colon word is k:g:i:j (g=graphlet Ordinal, i,j is a node (NOT ORBIT) pair in g.), followed by a count.
# We call (i,j) a "canonical node pair", or cnp for short.

# Bins need to be integers but we want reasonably high resolution, even when the sample sizes are small leading to
# low scores. The disadvantage of so many bins is that gawk takes up lots of RAM... oh well...
BINS_PER_WEIGHT_UNIT=100

hawk 'function WeightToBin(w){return int('$BINS_PER_WEIGHT_UNIT'*w);} # because weights are floats but we need to bin them
    BEGIN{windowSize='$WINDOW_SIZE'; '"$MINIMUMS"'; cNorm='$DO_NORM'} # cNorm = Boolean: normalize by StatMean(cnp)?
    ARGIND<=1+cNorm{
	uv=$1 # node pair
	ASSERT(2==split(uv,a,":"),"first column not colon-separated");
	u=a[1]; v=a[2];
	if(ARGIND==1) { # only need these assertions on first pass
	    ASSERT($2==0 || $2==1, "expecting second column to be Boolean");
	    #ASSERT(u>v,"u and v in wrong order");
	}
	if(u<v) {tmp=u; u=v; v=tmp}
	E[uv]=e[u][v]=$2 # edge Boolean
	for(i=3;i<NF;i+=2){ #col 3 onwards are (cnp,count) pairs
	    cnp=$i; c=$(i+1);
	    if(ARGIND==cNorm) StatAddSample(cnp,c); # if cNorm==0 this never gets executed
	    else {
		if(cNorm) c /= StatMean(cnp);
		#printf "Pearson[%s](%s,%g,%d)\n", uv, cnp,c,E[uv]
		PearsonAddSample(cnp,c,E[uv]);
	    }
	    c=WeightToBin(c); # integer because of histogram
	    #printf "bumping hist[%s][%g][%g]\n", cnp,c,E[uv]
	    ++hist[cnp][c][E[uv]]; # histogram: number of times orbit-pair $i had exactly count c for edge truth "e"
	    if(c>max[cnp])max[cnp]=c
	}
	#if(NR>100000)nextfile; #DEBUG
    }
    ENDFILE{printf "Finished reading ARGIND %d (%s)\n", ARGIND, FILENAME > "/dev/stderr"
	if(ARGIND==1+cNorm){
	    # Now produce empirical precision curve as a function of orbit-pair count (c), for each orbit-pair
	    printf("Predictive stats for '"$*"'\n") > "/dev/stderr"
	    for(cnp in hist) {
		keep=0 # is this cnp a keeper?
		if(_Pearson_N[cnp]>=min_samples && PearsonCompute(cnp) &&
		    ABS(_Pearson_rho[cnp])>=min_rho && _Pearson_t[cnp]>=min_t) {
		    #printf "cnp %s maxc = %d\n", cnp, max[cnp]
		    ASSERT(length(hist[cnp][max[cnp]])>0, "length(hist["cnp"][max_c="max[cnp]"])="length(hist[cnp][max[cnp]]));
		    numpr[cnp][ max[cnp] + 1 ]=0
		    for(c=max[cnp];c>=0; --c) { # starting at the highest orbit-pair counts
			numer[cnp]+=hist[cnp][c][1]; # those that actually had an edge
			denom[cnp]+=hist[cnp][c][1]+hist[cnp][c][0]; # both edge and non-edge
			# precision of this orbit-pair as func of c (add 1 to denom simply to avoid div by zero)
			if(denom[cnp]) prec[cnp][c] = numer[cnp]/denom[cnp];
			else prec[cnp][c] = prec[cnp][c+1];
			numpr[cnp][c]=numpr[cnp][c+1]+hist[cnp][c][0]; # number of genuine predictions from max_c down to c
			#printf "cnp %s c %d numer %d denom %d numpr %d prec %g\n", cnp, c, numer[cnp], denom[cnp], numpr[cnp][c], prec[cnp][c]
		    }
		    # Now, a heuristic "fix" to statistical noise for the few, highest-scoring counts: fix the top predictions
		    # so that the precision is artificially fixed to be monotonically increasing with orbit pair count
		    maxP[cnp]=0; # max "raw" precision
		    wPrec[cnp][0]=0;
		    for(c=0; c<=max[cnp];c++) { # starting at the LOWEST orbit-pair counts
			# make raw precision non-decreasing with increasing orbit count
			prec[cnp][c] = maxP[cnp] = MAX(maxP[cnp],prec[cnp][c]);
			# make weighted precision a moving average across windowSize values of c
			newPrec = (wPrec[cnp][c-1] * (windowSize-1) + prec[cnp][c])/windowSize; 
			wPrec[cnp][c] = MAX(newPrec, wPrec[cnp][c-1]); # precision should be monotonically increasing
			#printf "prec[%d] = %g; wPrec = %g; numpr = %d \n", c, prec[cnp][c], wPrec[cnp][c], numpr[cnp][c]
		    }
		    # re-purpose maxP to be the max wPrec (above it was max "raw" prec)
		    maxP[cnp]=wPrec[cnp][max[cnp]]; # maxP occurs at the highest c since wPrec is monotonically increasing
		    if(maxP[cnp] >= min_p) keep=1
		}

		if(keep) print PearsonPrint(cnp),cnp,maxP[cnp] > "/dev/stderr"
		else {
		    delete max[cnp]; delete maxP[cnp]; delete _Pearson_N[cnp];
		    delete numer[cnp]; delete denom[cnp]; delete prec[cnp]; delete wPrec[cnp]; delete hist[cnp];
		}
	    }
	}
	if('$PREDICTORS_ONLY')exit(0);
    }
    ARGIND==2+cNorm{ # same file, making predictions on the last pass
	uv=$1 # node pair
	ASSERT(2==split(uv,a,":"),"first column not colon-separated");
	u=a[1]; v=a[2]; if(u<v) {tmp=u; u=v; v=tmp}

	# First pass through the line: find the highest precision orbit pair
	c = p1 = 0;
	bestCol="none"; bestVal=0;
	for(i=3;i<NF;i+=2)if($i in hist) {
	    cnp=$i; c=$(i+1); if(cNorm) c /= StatMean(cnp)
	    c=WeightToBin(c);
	    if(c>max[cnp])c=max[cnp]; # clip the orbit count to the highest seen during training
	    # filter on "reasonable" orbit pairs
	    if(p1<wPrec[cnp][c]){p1=wPrec[cnp][c];bestCol=cnp;bestVal=c}
	}

	if(('$ALL_PREDICTIONS' || p1>min_p) && E[uv]<='$INCLUDE_KNOWN'){
	    printf "%s\t%g\t%d %s %d", uv, p1, bestVal, bestCol, numpr[bestCol][bestVal]
	    if('$VERBOSE') printf "\t[%s]\n",$0
	    else print ""
	}
    }' "$TMPDIR/input" "$TMPDIR/input" $THIRD_PASS |
    sort -k 2,2gr -k 3,3nr -k 5,5nr |
    if [ "$EVALUATE" = "" ]; then
	cat
    else
	awk 'ARGIND==1{E[$1][$2]=E[$2][$1]=1}
	    ARGIND==2{
		split($1,a,":"); T=((a[1] in E) && (a[2] in E[a[1]]));
		sum+=T; printf "%d %d prec %g %s\n", T, FNR, sum/FNR,$0
	    }' "$EVALUATE" -
    fi
