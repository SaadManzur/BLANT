#!/bin/sh
USAGE="$0 [-B blant_exe] [-s SampleMethod] [-M numSamples in millions] network.el [list of k values, defaults to just 4]
PURPOSE: Given nothing but an input graph, create edge preditions by running 'blant_exe -mp' and processing the output online.
(blant_exe defaults to './blant'.)
We cut the input network in half NODE-WISE (note that you lose 50% of all edges, but it's a more rigorous test than edgewise).
Then remove 10% of edges of both, giving two edge lists that you will attempt to predict. The result is 4 files:
Network G0.el plus G0-missing.el; and G1.el plus G1-missing.el.
Then run 'blant_exe -mp' on G0.el and G1.el, independently.
Then use the *histogram* from G0.el to applied to the blant-mp output of G1.el to predict G1-missing.el; then the reverse.
The idea is that we use G0.el to learn which motifs are predictive; then use that info on the blant-mp run of G1.el to
use the motif *info* in G1.el to rank predictions on missing edges from G1.el, then see how many of G1-missing.el we
recover."

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# generally useful Variables
NL='
'
TAB='	'

# Functions
die(){ (echo "USAGE: $USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ hawk "BEGIN{print $@}" </dev/null; }
integers() {
    case $# in
    1) awk 'BEGIN{for(i=0;i<'"$1"';i++) print i; exit}' ;;
    2) awk 'BEGIN{for(i='"$1"';i<='"$2"';i++) print i; exit}' ;;
    3) awk 'BEGIN{for(i='"$1"';i<='"$2"';i+='"$3"') print i; exit}' ;;
    esac
}
cpus() {
    TMP=/tmp/cpus.$$
    trap "/bin/rm -f $TMP; exit" 0 1 2 3 15

    # Most Linux machines:
    lscpu >$TMP 2>/dev/null && awk '/^CPU[(s)]*:/{cpus=$NF}END{if(cpus)print cpus; else exit 1}' $TMP && return

    # MacOS:
    ([ `arch` = Darwin -o `uname` = Darwin ] || uname -a | grep Darwin >/dev/null) && sysctl -n hw.ncpu && return

    # Cygwin:
    case `arch` in
    CYGWIN*) grep -c '^processor[ 	]*:' /proc/cpuinfo; return ;;
    *) if [ -d /dev/cpu -a ! -f /dev/cpu/microcode ]; then
	ls -F /dev/cpu | fgrep -c
	return
       fi
	;;
    esac

    # Oops
    echo "couldn't figure out number of CPUs" >&2; exit 1
}
randomizeLines() {
hawk 'BEGIN{Srand()}
    function randint(N){return int(N*rand())}
    {line[NR]=$0}
    END{
	N=NR;
	while(N) {
	    k=randint(N)+1;
	    print line[k];
	    line[k]=line[N--];
	}
    }' "$@"
}
induce() {
    FRAC=UNDEF
    INFILE=""
    [ -d $TMPDIR ] || die "induce expects TMPDIR '$TMPDIR' to already exist"
    mkdir -p $TMPDIR/induce
    case "$1" in
    -f) FRAC="$2"; shift 2;;
    *) [ -r "$1" ] || die "'$1' needs to be an input file of nodes";
	INFILE="$1"; shift;
	cp "$INFILE" $TMPDIR/induce/nodes-to-induce
	;;
    esac

    [ $# -eq 1 -a -r "$1" ] || die "Last argument must be a graph file"

    cp -p "$1" $TMPDIR/induce/graph.el || die "couldn't copy graph file '$1'"

    if [ "$FRAC" != "UNDEF" ]; then # create random induced subgraph
	awk 'BEGIN{srand(); n=0; FRAC='"$FRAC"';
	    if(FRAC < 0 || FRAC > 1) {
		printf "FRAC (%g) must be between 0 and 1\n", FRAC >"/dev/fd/2";
		exit 1
	    }
	}
	{
	    if(!($1 in V)){ nodeName[++n]=$1; V[$1]=1}
	    if(!($2 in V)){ nodeName[++n]=$2; V[$2]=1}
	}
	END {
	    n1=int(FRAC*n);
	    induced[""]=1; delete induced[""]; # create empty array
	    for(i=0;i<n1;i++) {
		do {u = nodeName[int(1+n*rand())]} while(u in induced);
		induced[u]=1;
		print u;
	    }
	}' $TMPDIR/induce/graph.el > $TMPDIR/induce/nodes-to-induce
    fi
    [ $? -eq 0 ] || die "FRAC ($FRAC) not valid"

    #ls -l $TMPDIR/induce/nodes-to-induce $TMPDIR/induce/graph.el >&2
    awk 'ARGIND==1{
	    if(NF!=1){
		printf "first file \"%s\", line %d: expecting 1 column (a list of nodes) but got \"%s\"\n",
		    FILENAME,FNR,$0 >"/dev/fd/2"; exit 1
	    }
	    nodes[$1]=1
	}
	ARGIND==2{
	    if(NF!=2) {
		printf "second file \"%s\", line %d: expecting 2 columns (an edge) but got \"%s\"\n",
		    FILENAME,FNR,$0 >"/dev/fd/2"; exit 1
	    }
	    if(($1 in nodes) && ($2 in nodes))print;
	}' $TMPDIR/induce/nodes-to-induce $TMPDIR/induce/graph.el
}

##### Now the actual code #######
S=MCMC
N=1
Z=000000 # default is millions (6 zeros)
BLANT=./blant
while true; do
    case "$1" in
    -[Zz]) Z="$2"; shift 2;;
    -[Mm]) N="$2"; shift 2;;
    -[Ss]) S="$2"; shift 2;;
    -[Bb]) BLANT="$2"; shift 2;;
    -[Tt][Mm][Pp]*) TMPDIR="$2"; shift 2;;
    -) break;; # using stdin as input
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

[ "$1" = - -o -s "$1" ] || die "expecting an edgeList.el file as input (stdin is allowed but must be specified as '-')"
echo "Input file is '$1'" >&2

[ "$TMPDIR" = "" ] && TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
echo "TMPDIR is $TMPDIR" >&2

[ -s $TMPDIR/G.el ] || cat "$1" | tee $TMPDIR/G.el | newlines | sort -u | tee $TMPDIR/nodes | randomizeLines | awk '{m=(NR-1)%2; print > "'$TMPDIR'/G"m".nodes"}'

shift # get rid of the filename argument

for i in 0 1; do
    [ -s "$TMPDIR/G$i.el" ] || induce $TMPDIR/G$i.nodes $TMPDIR/G.el > $TMPDIR/G$i.el
    m=`wc -l <$TMPDIR/G$i.el | awk '{print int($1/10)}'`
    [ -s $TMPDIR/G$i-train.el ] || head -$m $TMPDIR/G$i.el | tee $TMPDIR/G$i-test.el | fgrep -v -f - $TMPDIR/G$i.el > $TMPDIR/G$i-train.el
done

CORES=`cpus`
echo "Using at most $CORES cores (or set PARALLEL to bash if you have trouble with parallel)"
PARALLEL="parallel $CORES"
# PARALLEL=sh # use this if parallel doesn't work.
[ $# -eq 0 ] && set 4
echo "Using k values" "$@"
for k
do
    [ 4 -le $k -a $k -le 8 ] || die "k values must be between 4 and 8"
    threads=1
    maxExpense=100000000 # 100 million
    expense=`parse "$N$Z*fact($k)/fact(4)"` # approximate expense as a function of k
    [ "$expense" -gt $maxExpense ] && threads=`parse "int($expense/$maxExpense)"`
    samplesPerThread=`parse "int($N$Z/$threads)"`
    for DUP in `integers $threads`; do
	for i in 0 1; do
	    OUT="$TMPDIR/blant-mp.train$i-$N$Z.$S.k$k.$DUP"
	    [ -s "$OUT" ] || echo "time $BLANT -mp -s $S -n $samplesPerThread -k $k $TMPDIR/G$i-train.el | tee '$OUT' | wc"
	done
    done
done > $TMPDIR/jobs.txt
[ -f $TMPDIR/ABORT ] && die "`cat $TMPDIR/ABORT`"
echo "Starting `wc -l < $TMPDIR/jobs.txt` jobs from $TMPDIR/jobs.txt running at most $CORES in parallel" >&2
cat $TMPDIR/jobs.txt | $PARALLEL
for i in 0 1; do
    OUT="$TMPDIR/blant-mp.train$i-$N$Z.$S.K"
    if [ -s "$OUT" ]; then
	echo "$OUT exists and has nonzero size; skipping" >&2
    else
	if [ `cd $TMPDIR; ls blant-mp.train$i-*.k* | wc -l` -eq 1 ]; then
	    (cd $TMPDIR; ln -s blant-mp.train$i-*.k* blant-mp.train$i-$N$Z.$S.K)
	else
	    (cd $TMPDIR; echo Merging blant-mp.train$i-*.k*) >&2
	    echo "./predict-merge.sh $TMPDIR/blant-mp.train$i-*.k* > $TMPDIR/blant-mp.train$i-$N$Z.$S.K"
	fi
    fi
done | $PARALLEL

for i in 0 1; do
    ./predict-blant-mp.sh "$TMPDIR/blant-mp.train$i-$N$Z.$S.K" | cut -f1 | sed 's/:/	/' | #hawk '{gsub(":"," "); printf "%s\t%s\n",MAX($2,$3),MIN($2,$3)}' |
    fgrep -f - $TMPDIR/G$i-test.el | 
    awk 'END{m=2*'$m'; printf "%d correct predictions out of %d, recall = %g\n",NR,m,NR/m}'
