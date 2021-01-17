#!/bin/sh
USAGE="$0 [-B blant_exe] [-s SampleMethod] [-M numSamples in millions (default 10)] [-include-known] network.el [k values, default 4 only]
PURPOSE: Given nothing but an input graph, create edge preditions by running 'blant_exe -mp' and processing the output online.
(blant_exe defaults to './blant'.). By default the predictions only include edges not in the original network (ie.,
genuine predictions), but if '-include-known' is specified, the 'predictions' include edges from the original network
as well, which facilitates evaluating the predictions on the training set."

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# generally useful Variables
NL='
'
TAB='	'

# Bourne shell Functions
die(){ (echo "USAGE: $USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ hawk "BEGIN{print $@}" </dev/null; }
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
    die "couldn't figure out number of CPUs"
}

##### Now the actual code #######
S=MCMC
M=1
Z=000000 # number of zeroes after a million
BLANT=./blant
INCLUDE_KNOWN=''
while true; do
    case "$1" in
    -[Mm]) M="$2"; shift 2;;
    -[Ss]) S="$2"; shift 2;;
    -[Bb]) BLANT="$2"; shift 2;;
    -[Tt][Mm][Pp]*) TMPDIR="$2"; shift 2;;
    -include-known) INCLUDE_KNOWN="$1"; shift;;
    -) break;; # using stdin as input
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

[ "$TMPDIR" = "" ] && TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
export TMPDIR
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

[ "$1" = - -o -s "$1" ] || die "expecting an edgeList.el file as input (stdin is allowed but must be specified as '-')"
echo "Input file is '$1'" >&2
cp -p "$1" $TMPDIR/G.el
shift # get rid of the filename argument
echo "TMPDIR is $TMPDIR" >&2

CORES=`cpus`
#echo "Using at most $CORES cores (or set PARALLEL to bash if you have trouble with parallel)" >&2
PARALLEL="parallel $CORES"
# PARALLEL=sh # use this if parallel doesn't work.
[ $# -eq 0 ] && set 4
#echo "running blant on the following jobs using k values" "$@" >&2
for k
do
    [ 4 -le $k -a $k -le 8 ] || die "k values must be between 4 and 8"
    threads=1 # below is the approximate relative expense compared to k=4
    #expense=`echo $k | awk '{prod=1; K=$1; for(k=4;k<K;k++)prod *= sqrt(2)*(k-1)*(k-2)*(k-3); print int(prod)}'`
    #[ "$expense" -gt $M ] && threads=`parse "MAX(1,int($expense/$M/100))"`
    samplesPerThread=`parse "int($M$Z/$threads)"`
    #echo "Relative expense is $expense with so $threads threads needed for k=$k" >&2
    for DUP in `seq 1 $threads`; do
	OUT="$TMPDIR/blant-mp.$S.${M}M.k$k.$DUP"
	[ -s "$OUT" ] || echo "time $BLANT -mp -s $S -n $samplesPerThread -k $k $TMPDIR/G.el | tee '$OUT' | wc >&2"
    done
done > $TMPDIR/jobs.txt
[ -f $TMPDIR/ABORT ] && die "`cat $TMPDIR/ABORT`"
echo "Starting `wc -l < $TMPDIR/jobs.txt` jobs from $TMPDIR/jobs.txt running at most $CORES in parallel using k values $@" >&2
cat $TMPDIR/jobs.txt | $PARALLEL
OUT="$TMPDIR/blant-mp.$S.${M}M.K"
if [ -s "$OUT" ]; then
    echo "$OUT exists and has nonzero size; skipping" >&2
else
    if [ `ls   $TMPDIR/blant-mp.*.k?.* | wc -l` -eq 1 ]; then
	(cd $TMPDIR; echo Only one file, linking blant-mp.*.k?.* "->" blant-mp.$S.${M}M.K >&2
	ln -s                       blant-mp.*.k?.*      blant-mp.$S.${M}M.K)
    else
	(cd $TMPDIR; echo Merging blant-mp.*.k?.*) >&2
	#./predict-merge.sh $TMPDIR/blant-mp.*.k?.*                               > $TMPDIR/blant-mp.$S.${M}M.K
	 cat                $TMPDIR/blant-mp.*.k?.* | $BLANT -k8 -mq $TMPDIR/G.el > $TMPDIR/blant-mp.$S.${M}M.K
    fi
fi
./predict-from-blant-mp.sh $INCLUDE_KNOWN "$TMPDIR/blant-mp.$S.${M}M.K" | cut -f1
