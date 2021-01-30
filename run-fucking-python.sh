#/bin/sh
# Do whatever is necessary to load your favorite version of python on the line(s) below
module load python/3.6.8

# Now tell me the name of the python interpreter that you want to use, given the above
PYTHON=python3

# PYTHON GENIUSES DO NOT TOUCH ANYTHING BETWEEN HERE, AND WHERE IT SAYS "Put your fucking Python shit ..."
#############################################################################################################

USAGE="$0 bla bla bla
PURPOSE: Run a python script included at the end of the Bourne Shell script."

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

sed '1,/^#Put/d' "$0" > $TMPDIR/fuck-python.py
chmod +x $TMPDIR/fuck-python.py

exec $PYTHON $TMPDIR/fuck-python.py "$@"
exit 1


#############################################################################################################
#Put your fucking Python shit below this line.

print("Hello World")
