#include "Oalloc.h" // "Once alloc", for allocating tons of little integers via pointers
#include "../blant.h"
#include "../blant-output.h"
#include "../blant-utils.h"
#include "blant-predict.h"
#include "htree.h"
#include "rand48.h"
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include <unistd.h>

// Watch memory usage
#include <sys/time.h>
#include <sys/resource.h>

// Reasonable Defaults for a Mac with 16GB RAM since it's too much trouble to find the real values.
static double MAX_GB, totalGB = 16; //, freeGB = 8, MAX_GB;

#define GB (1024.0*1024.0*1024.0)
typedef sig_t sighandler_t;
#if __APPLE__
  #define RUSAGE_MEM_UNIT 1L
  #define SYSINFO_MEM_UNIT 1L
#else
  #include <sys/sysinfo.h>
  #if __WIN32__ || __CYGWIN__
    #define SYSINFO_MEM_UNIT 4096L  // 4k unit???
    #define RUSAGE_MEM_UNIT 4096L  // 4k page
  #else
    #define SYSINFO_MEM_UNIT info.mem_unit
    #define RUSAGE_MEM_UNIT info.mem_unit
  #endif
#endif

void CheckRAMusage(void)
{
    static struct rusage usage;
    static int status;
#if !__APPLE__
    static struct sysinfo info;
    status = sysinfo(&info);
    assert(status == 0);
    //Note("RAW sysinfo data: totalram %ld freeram %ld\n", info.totalram, info.freeram);
    totalGB = info.totalram*1.0*SYSINFO_MEM_UNIT/GB; 
    //freeGB = info.freeram*1.0*SYSINFO_MEM_UNIT/GB;
#endif
    MAX_GB=totalGB-4; //MIN(totalGB/2, freeGB-4); // at most half the machine's RAM, and leaving at least 4GB free
    if(!usage.ru_maxrss) Note("System claims to have totalram %g GB;  aiming to use MAX %g GB",
	totalGB, MAX_GB);
    status = getrusage(RUSAGE_SELF, &usage);
    assert(status == 0);
    //Note("RAW rusage data: res %ld drss %ld srss %ld", usage.ru_maxrss, usage.ru_idrss, usage.ru_isrss);
    if(usage.ru_maxrss*RUSAGE_MEM_UNIT/GB > MAX_GB || (usage.ru_idrss+usage.ru_isrss)*RUSAGE_MEM_UNIT/GB > MAX_GB)
    {
	double new=usage.ru_maxrss*(double)RUSAGE_MEM_UNIT/GB;
	static double previous;
	if(new > 1.1*previous) {
#if !__APPLE__
	    status = sysinfo(&info);
	    //freeGB = info.freeram*SYSINFO_MEM_UNIT/GB;
#endif
	    Warning("WARNING: Resident memory usage has reached %g GB", new);
	    previous = new;
	}
	_memUsageAlarm=true;
    }
}

/* General note: when associating any pair of nodes (u,v) in the large graph G, with "canonical orbit pairs" in
** canonical graphlets, the number of possible combination is so incredibly enormous and sparse, that even hashing
** was too wasteful of memory. Thus I moved to storing the associations for *each* pair (u,v) in one binary tree.
** Thus for G with n nodes, there are potentially (n choose 2) distinct binary trees, each holding the list
** of canonical orbit pair associations and the respective count. Basically for a 10,000 node graph, and using
** k=8, there are almost a quarter million  possible orbit pairs... so the total number of possible *associations*
** between any (u,v) in G and any canonical node pair is (n choose 2) * 244,000 ... and for n=10,000 nodes, that's
** 1.2 x 10^13 associations... each of which must have a *count* (potentially requiring 32 bits), for a total
** of 5e13 bytes=50 TB. So not feasible, and even a hash table (which is sparsely stored to reduce collisions)
** was requiring >100GB to store these associations. A Binary tree is a bit slower but requires far less RAM.
*/

// Note for the PredictGraph we actually use an HTREE--heirarchical Tree, where the first key is the
// orbit pair (o,p) from graphlet K, and the second key is any edge in any graphlet K; we then use the
// total number of edges from G seen in any graphlet K to compute the weight of u:v's o:p orbit pair.
static HTREE ***_PredictGraph; // lower-triangular matrix (ie., j<i, not i<j) of dictionary entries
static BINTREE *_predictiveOrbits; // if non-NULL, any orbit not in this dictionary is ignored.

// Is is there anything in _PredictGraph[i][j]?
#define PREDICT_GRAPH_NON_EMPTY(i,j) (_PredictGraph[i][j] && _PredictGraph[i][j]->tree->root)

// This is a hack that's required for parallelism: technically we only need to output the *final* counts once
// we've accumulated them. But if we're a child process just shunting those counts off to a parent process who's
// accumulating these counts globally, it turns out that it's expensive for the parent, and if we wait until the
// end to provide our output then the parent is sitting around twiddling its thumbs until the end, and is suddenly
// inundated with expensive parsing. So instead, if we're a child process, every fraction of a second (0.1s seems
// best), we spit out our accumulation so far so the parent can parse them online. This Boolean is set to true each
// 0.1s, and below we check it to see if it's time to flush our counts.
static Boolean _flushCounts = true;

// Signal handler for the SIGALRM that occurs periodically forcing us to flush our counts to a parent (if we're _child).
int AlarmHandler(int sig)
{
    //assert(_child && !_flushCounts);
    //fprintf(stderr, "Alarm in process %d\n", getpid());
    alarm(0);
    CheckRAMusage();
    if(_child) _flushCounts = true;
    signal(SIGALRM, (sighandler_t) AlarmHandler);
    alarm(1); // set alarm for 1 second hence
    return 0;
}

// Allocate the NULL pointers for just the *rows* of the PredictGraph[i].
void Predict_Init(GRAPH *G) {
    int i;
    assert(_PredictGraph==NULL);
    // This just allocates the NULL pointers, no actual binary trees are created yet.
    _PredictGraph = Calloc(G->n, sizeof(HTREE**)); // we won't use element 0 but still need n of them
    for(i=1; i<G->n; i++) _PredictGraph[i] = Calloc(i, sizeof(HTREE*));

    char *predictive = getenv("PREDICTIVE_ORBITS");
    if(predictive) {
	int numPredictive=0;
	_predictiveOrbits = BinTreeAlloc((pCmpFcn) strcmp, (pFointCopyFcn) strdup, (pFointFreeFcn) free, NULL, NULL);
	// Is it a string containing orbits, or a filename?
	char *s = predictive;
	// The only characters allowed in orbit lists are digits, spaces, and colons.
	while(*s && (isdigit(*s) || isspace(*s) || *s == ':')) s++;
	if(*s) { // we found an invalid orbit-descriptor character; assume it's a filename
	    FILE *fp = Fopen(predictive,"r");
	    char buf[BUFSIZ];
	    while(1==fscanf(fp,"%s",buf)) {
		++numPredictive;
		BinTreeInsert(_predictiveOrbits, (foint) buf, (foint) NULL);
	    }
	    fprintf(stderr, "Read %d predictive orbits from file %s\n", numPredictive, predictive);
	} else { // we got to the end of the string, it's a raw list of orbits
	    FILE *fp = fopen(predictive,"r");
	    if(fp) Fatal("PREDICTIVE_ORBITS <%s> looks like a list of orbits but there's also a file by that name", predictive);
	    fprintf(stderr, "Reading $PREDICTIVE_ORBITS: <%s>\n", predictive);
	    s=predictive;
	    while(*s)
	    {
		while(isspace(*s)) s++; // get past whitespace
		if(!*s) break;
		predictive = s;
		while(*s && !isspace(*s)) s++; // get to next whitespace
		if(isspace(*s)) *s++ = '\0'; // terminate the string
		++numPredictive;
		BinTreeInsert(_predictiveOrbits, (foint) predictive, (foint) NULL);
	    }
	    fprintf(stderr, "Read %d orbits from $PREDICTIVE_ORBITS\n", numPredictive);
	}
    }

    if(_child) {
	signal(SIGALRM, (sighandler_t) AlarmHandler);
	alarm(1);  // (int)(2*_MAX_THREADS*drand48())); // on average, have one child per second send data to the parent process
    }
}


/*
** Our ultimate goal is: for a pair of nodes (u,v) and an orbit pair (o,p), count the number of distinct quadruples
** (q,r,x,y) where (q,r)!=(o,p) is an edge, and (x,y) is an edge in G also !=(u,v). NOTE: after some experimentation
** it may be sufficient to ignore (q,r) and only count the edges (x,y) in G.
**
** In the case that we can ignore (q,r), it's fairly easy: our "count only" version creates one sorted binary
** tree for each (u,v) pair in G, with the tree sorted on the key (o,p) represented as a string, and with the data
** member being an integer count. Adding (x,y) into the mix is easy: instead of the data member being a count, it'll
** be *another* binary tree using the key (x,y), and *that* binary tree will have no data member. Then the count for
** (u,v,o,p) is simply the number of entries in the (u,v)->(o,p) binary tree, ie the number of unique (x,y) keys in it.
**
** To use a binary tree for the more complex case, we *could* do it as follows: for each (u,v) and each (o,p), have a
** separate binary tree--which multiplies our number of binary tree by (k choose 2)--and then in each such binary tree,
** keep track of the number of distinct keys (q,r,x,y); we wouldn't even need a data member, we just need to count
** distinct keys, which would be B->n. [THIS IS HOW WE ACTUALLY DO IT, using an HTREE.]

** If we wanted to be more clever, we could do some clever encoding. (We won't unless saving memory becomes crucial.)
** Note that (o,p,q,r) has exactly (k choose 2)^2 possible values [or just (k choose 2) if we ignore (q,r).].
** At k=8 that's only 784 [28] possible values, and we only need to remember a BOOLEAN of each.
** Thus, we'll represent whether we've seen (o,p,q,r) as a SET*, which will require
** about 100 [4] bytes total. Given (o,p,q,r) we'll convert (o,p) to int via creating an empty TinyGraph, adding (o,p),
** then op=TinyGraphToInt; same with (q,r) giving qr; finally opqr=op*(k choose2) + qr.
** Then, to fully encode the (u,v,o,p,q,r,x,y) octuplet [sextuplet], we'll keep the _PredictGraph[u][v]
** binary trees, but now the *key* will simply be "x:y", the "internal edge", and then the octuplet [sextuplet]
** can be queried as: key = BinaryTreeLookup(PG[u][v], "x:y"); SET *uvxy=(SET*) key; and finally SetAdd(uvxy,opqr).
** Finally retrieving and printing the output will require, for each (u,v) pair, to traverse its binary
** tree across all (x,y) pairs, accumulating sum[qr] += !!SetIn(uvxy, opqr) (!! to ensure it's 0 or 1).
** In English, that's saying: for a given (u,v) pair, its value at CNP (o,p) is the sum, across all (x,y)
** edges that have been observed in the same sampled graphlet, of [1 if ignoring qr, or] whether (x,y) has
** appeared at CNP (q,r). (1 if yes, 0 if no).
*/
#define COUNT_uv_only 1 // this works OK but COUNT_xy_only works better
#if COUNT_uv_only
  #define DEG_WEIGHTED_COUNTS 1
#else
  #define COUNT_xy_only 1 // count unique [xy] edges only; othewise include both q:r and x:y; edges only seems to work best.
#endif

#ifdef COUNT_xy_only
// This function is only used if we want to see how many times each xy edge has been seen--but it doesn't seem necessary.
Boolean Traverse_xy(foint key, foint data) {
    if(!COUNT_xy_only) Apology("!COUNT_xy_only not yet implemented");
    char *ID = key.v;
    printf("(x:y)=%s %g",ID, data.f);
    return true;
}
#endif

/* TraverseNodePairCounts is called indirectly via BinTreeTraverse, which prints *all* the currently stored
** participation counts of a particular pair of nodes (u,v) in G. This function is called on each such
** participation count, and it simply prints the k:o:p canonical orbit pair, and its weight (|xy|), to stdout.
*/
Boolean TraverseNodePairCounts(foint key, foint data) {
    char *ID = key.v;
#if COUNT_uv_only
    printf("\t%s %g",ID, data.f);
#else // seems to be the same output whether or not COUNT_xy_only, since it's the ID that distinguishes them.
    BINTREE *op_xy = data.v;
    printf("\t%s %d",ID, op_xy->n);
    //BinTreeTraverse(op_xy, Traverse_xy); // only if we want to the count of each edge or (q,r,x,y) quad.
#endif
    return true;
}


/* Loop across all pairs of nodes (u,v) in G, and for each pair, print one line that contains:
** the pair of nodes, the edge truth, and all the participation counts.
*/
void PredictFlushAllCounts(GRAPH *G){
    alarm(0); // turn off _flushCounts alarm
    int i,j;
    for(i=1; i < G->n; i++) for(j=0; j<i; j++) {
	if(PREDICT_GRAPH_NON_EMPTY(i,j))  // only output node pairs with non-zero counts
	{
	    printf("%s %d", PrintNodePairSorted(i,':',j), GraphAreConnected(G,i,j));
	    BinTreeTraverse(_PredictGraph[i][j]->tree, TraverseNodePairCounts);
	    puts("");
	}
    }
}


static void UpdateNodePair(int G_u, int G_v, char *ID, double count); // Prototype; code is below


// TraverseCanonicalNodePairs is called via BinaryTreeTraverse to traverse all the known associations between a particular pair
// of nodes in a canonical graphlet, in order to transfer that information to a a pair of nodes from a sampled graphlet in G.
// (A new traversal is constructed for every pair of nodes in the sampled graphlet.)
// This functions' job is to transfer that info to global node pairs in G called (G_u,G_v).
// The pair of nodes (G_u,G_v) in G must be global since there's no other mechanism to get that info here.
static int _TraverseCanonicalPairs_G_u, _TraverseCanonicalPairs_G_v; // indices into _PredictGraph[G_u][G_v]
static char *_TraverseCanonicalPairs_perm; // canonical permutation array for the current sampled graphlet...
static double _TraverseCanonicalPairs_deg_weight; // degree-based weigt of interior nodes
static unsigned *_TraverseCanonicalPairs_Varray; // ...and the associated Varray.

static Boolean TraverseCanonicalPairs(foint key, foint data) {
    char *ID = key.v; // ID is a *string* of the form k:o:p[[q:r]:x:y], where (o,p) is the canonical orbit pair.
    int *pCanonicalCount = data.v; // the count that this (o,p) pair was seen at (u,v) during canonical motif enumeration.
    static char ID2[BUFSIZ], ID3[BUFSIZ];
#if COUNT_uv_only // the tail end of the ID string contains either x:y (COUNT_xy_only) or q:r:x:y (otherwise)
    strcpy(ID3, ID); // make a copy so we can nuke bits of it.
#else // the tail end of the ID string contains either x:y (COUNT_xy_only) or q:r:x:y (otherwise)
    strcpy(ID2, ID); // make a copy so we can nuke bits of it.
    int q,r, IDlen = strlen(ID);
    q = *(ID+IDlen-3)-'0';
    r = *(ID+IDlen-1)-'0';
    assert(0<=q && 0 <= r && q<_k && r<_k);
    // x and y in the non-canonical motif g that is induced from G
    int g_x=_TraverseCanonicalPairs_perm[q], g_y=_TraverseCanonicalPairs_perm[r];
    int G_x=_TraverseCanonicalPairs_Varray[g_x], G_y=_TraverseCanonicalPairs_Varray[g_y]; // x and y in BIG graph G.
    if(g_x < g_y) { int tmp = g_x; g_x=g_y; g_y=tmp; } // lower triangle of g
    if(G_x < G_y) { int tmp = G_x; G_x=G_y; G_y=tmp; } // for consistency in accessing _PredictGraph
    char *pColon = (ID2+IDlen-4); // prepare to nuke the : before q:r
    assert(*pColon == ':'); *pColon = '\0';
  #if COUNT_xy_only
    sprintf(ID3, "%s:%d:%d", ID2, G_x,G_y); // becomes k:o:p:x:y, where x and y are from G (not g)
  #else
    sprintf(ID3, "%s:%d:%d:%d:%d", ID2, q,r,G_x,G_y); // k:o:p:q:r:x:y
  #endif
#endif
    UpdateNodePair(_TraverseCanonicalPairs_G_u, _TraverseCanonicalPairs_G_v, ID3,
	*pCanonicalCount * _TraverseCanonicalPairs_deg_weight);
    return true;
}


// Given a pair of nodes (u,v) in G and an association ID, increment the (u,v) count of that association by the count
// of the ID. Note that unless COUNT_uv_only is true, this function is *only* used during merge mode (-mq).
static void UpdateNodePair(int G_u, int G_v, char *ID, double count) {
    foint fUVassoc; // either a count, or a pointer to either a count or sub-binary tree.
    if(G_u<G_v) { int tmp=G_u; G_u=G_v;G_v=tmp;}
    if(_PredictGraph[G_u][G_v] == NULL)
	_PredictGraph[G_u][G_v] = HTreeAlloc(1+!COUNT_uv_only); // depth 1 if COUNT_uv_only, or 2 for any COUNT_*xy
    char buf[BUFSIZ], *IDarray[2];
    IDarray[0] = strcpy(buf,ID);
#if !COUNT_uv_only
    // find the 3rd colon so we can isolate k:o:p from [q:r:]x:y
    int i, G_x,G_y;
    char *s=buf;
    for(i=0;i<3;i++)while(*s++!=':') ;
    assert(*--s==':'); *s++='\0';
#if PARANOID_ASSERTS && COUNT_xy_only
    assert(2==sscanf(s,"%d:%d",&G_x,&G_y)); assert(G_x > G_y);
#endif
    IDarray[1] = s; // the string [q:r:]x:y
#endif
    if(HTreeLookup(_PredictGraph[G_u][G_v], IDarray, &fUVassoc))
	fUVassoc.f += count;
    else 
	fUVassoc.f = count;
    HTreeInsert(_PredictGraph[G_u][G_v], IDarray, fUVassoc);
#if PARANOID_ASSERTS
    float newValue = fUVassoc.f;
    assert(HTreeLookup(_PredictGraph[G_u][G_v], IDarray, &fUVassoc)
	&& fabs(fUVassoc.f - newValue)/newValue < 1e-6);
#endif
    //printf("P %s %s %g (%g)\n", PrintNodePairSorted(G_u,':',G_v), ID, count, fUVassoc.f);
}


/* This function is used to merge the types of lines above (u:v edge, k:g:i:j count, etc) from several
** sources. It processes one line to figure out the node pair, and all the participation counts, and 
** then adds that data to the internal tally, which will later be printed out. Used either when -mp mode
** is invoked with multi-threading (-t N) or in predict_merge mode (-mq).
** NOTE ON EFFICIENCY: it may be more efficient to use built-in functions like strtok or scanf but I
** couldn't esaily get them to work so decided to hack it together myself. It's NOT particularly efficient
** and could probably be improved. (I think the BinaryTree stuff is not the bottleneck, it's the code
** right here.)
*/
void Predict_ProcessLine(GRAPH *G, char line[])
{
    assert(!_child);
    if(line[strlen(line)-1] != '\n') {
	assert(line[strlen(line)-1] == '\0');
	Fatal("char line[%s] buffer not long enough while reading child line in -mp mode",line);
    }
    char *s0=line, *s1=s0;
    foint fu, fv;
    int G_u, G_v;
    char *nameGu, *nameGv;
    if(_supportNodeNames) {
	while(*s1!=':') s1++; *s1++='\0'; nameGu=s0; s0=s1;
	while(*s1!=' ') s1++; *s1++='\0'; nameGv=s0; s0=s1;
	if(!BinTreeLookup(G->nameDict, (foint)nameGu, &fu)) Fatal("PredictMerge: node name <%s> not in G", nameGu);
	if(!BinTreeLookup(G->nameDict, (foint)nameGv, &fv)) Fatal("PredictMerge: node name <%s> not in G", nameGv);
	G_u = fu.i; G_v = fv.i;
	//printf("Found names <%s> (int %d) and <%s> (int %d)\n", nameGu, G_u, nameGv, G_v);
    }
    else {
	while(isdigit(*s1)) s1++; assert(*s1==':'); *s1++='\0'; G_u=atoi(s0); s0=s1;
	while(isdigit(*s1)) s1++; assert(*s1==' '); *s1++='\0'; G_v=atoi(s0); s0=s1;
    }
    assert(0 <= G_u && G_u < G->n);
    assert(0 <= G_v && G_v < G->n);
    assert(*s0=='0' || *s0=='1');  // verify the edge value is 0 or 1
    // Now start reading through the participation counts. Note that the child processes will be
    // using the internal INTEGER node numbering since that was created before the ForkBlant().
    s1=++s0; // get us to the tab
    while(*s0 == '\t') {
	int kk,g,i,j, count;
	char *ID=++s0; // point at the k value
	assert(isdigit(*s0)); kk=(*s0-'0'); assert(3 <= kk && kk <= 8); s0++; assert(*s0==':');
	s0++; s1=s0; while(isdigit(*s1)) s1++; assert(*s1==':'); *s1='\0'; g=atoi(s0); *s1=':'; s0=s1;
	s0++; assert(isdigit(*s0)); i=(*s0-'0'); assert(0 <= i && i < kk); s0++; assert(*s0==':');
	s0++; assert(isdigit(*s0)); j=(*s0-'0'); assert(0 <= j && j < kk); s0++; assert(*s0==' '); *s0='\0';
	s0++; s1=s0; while(isdigit(*s1)) s1++;
	if(!(*s1=='\t' || *s1 == '\n'))
	    Fatal("(*s1=='\\t' || *s1 == '\\n'), line is \n%s", line);
	// temporarily nuke the tab or newline, and restore it after (need for the top of the while loop)
	char tmp = *s1;
	*s1='\0'; count=atoi(s0);
	assert(0 <= g && g < _numCanon && 0<=i&&i<kk && 0<=j&&j<kk);
	UpdateNodePair(MAX(G_u,G_v) , MIN(G_u,G_v), ID, count);
	*s1 = tmp;
	assert(*(s0=s1)=='\n' || *s0 == '\t');
    }
    assert(*s0 == '\n');
}


/* Rather than explicitly enumerating the sub-motifs of every graphlet sampled from G, we do the recursive enumeration
** of all motifs under a graphlet only once---on the canonical graphlet. Then we *store* the association between the
** nodes of the canonical graphlet, and the motif associations. The associations get stored in a binary tree, one
** for each pair of canonical nodes in each canonical graphlet. Then, when we sample a graphlet from G, we simply
** determine what canonical it is using our standard BLANT lookup table _K, and use the perm[] matrix to transfer
** the canonical associations to the nodes of the graphlet sampled from G. The information is transferred from
** canonical to Varray in the function "TraversCanonicalNodes" elsewhere in the code.
**
** For each canonical graphlet, and for each pair of nodes in the graphlet, we maintain a dictionary (as a binary tree)
** of all the participation counts for all motif node-pairs that that canonical graphlet node pair participates in,
** across all sub-motifs of that canonical graphlet.
** Note that WE DO NOT KNOW HOW TO PROPERTLY NORMALIZE YET.... so for now do the easiest thing, which is count them all.
*/
static BINTREE *_canonicalParticipationCounts[MAX_CANONICALS][MAX_K][MAX_K];

/*
** All we need to add here is the inner loop that you have further down, over q and r, storing
** the count in relation to o and p. Then, when you TraverseCanonicals later, you'll use _perm[]
** to recover both u:v and x:y, and then do one of the following:
**      if we want the weight of u:v:o:p to be unique q:r:x:y quadruples, then:
	    
** contain the count of unique edges that contribute to u:v.
*/
void SubmotifIncrementCanonicalPairCounts(int topOrdinal, int Gint, TINY_GRAPH *g)
{
#if PARANOID_ASSERTS
    assert(TinyGraph2Int(g,_k) == Gint);
#endif
    int i, j, GintOrdinal=_K[Gint];
    char perm[MAX_K];
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);
    for(i=1;i<_k;i++) for(j=0;j<i;j++) // all pairs of canonical nodes
    {
	// We are trying to determine the frequency that a pair of nodes in the topOrdinal have an edge based on their
	// being located at a pair of canonical nodes in a sub-motif. The frequency only makes sense if the underlying
	// edge between them can sometimes exist and sometimes not; but if TinyGraphAreConnected(u,v)=true, the edge
	// already exists in this motif (and thus also in the topOrdinal) and there's nothing to predict. Thus, we only
	// want to check this node pair if the motif does NOT have the edge. (Comment label: (+++))
	int u=(int)perm[i], v=(int)perm[j]; // u and v in the *current* motif
        if(!TinyGraphAreConnected(g,u,v))
	{
	    int o=_orbitList[GintOrdinal][i], p=_orbitList[GintOrdinal][j];
	    // the association is between a node *pair* in the canonical top graphlet, and a node *pair* of the
	    // motif that they are participating in. Since the pair is undirected, we need to choose a unique order
	    // to use an a key, so we simply sort the pair.
	    // And since we want lower triangle, row > column, always, so o=1...n-1, p=0...o-1
	    if(u < v) { int tmp = u; u=v; v=tmp; }
	    if(o < p) { int tmp = o; o=p; p=tmp; }
	    assert(_canonicalParticipationCounts[topOrdinal][u][v]);
	    int *pcount;
	    char buf[BUFSIZ], buf_kop_only[BUFSIZ];
	    sprintf(buf_kop_only,"%d:%d:%d", _k,o,p);
#if COUNT_uv_only
	    sprintf(buf,"%d:%d:%d", _k,o,p);
#else
	    int q,r;
	    for(q=1;q<_k;q++) for(r=0;r<q;r++) if(q!=o || r!=p) {
		int x=(int)perm[q], y=(int)perm[r];
		if(TinyGraphAreConnected(g,x,y)) { // (x,y) only counts if it's an edge
  #if COUNT_xy_only
		    sprintf(buf,"%d:%d:%d:%d:%d", _k,o,p,x,y);
  #else
		    sprintf(buf,"%d:%d:%d:%d:%d:%d:%d", _k,o,p,q,r,x,y);
  #endif
		}
	    }
#endif
	    if(!_predictiveOrbits || BinTreeLookup(_predictiveOrbits, (foint)(void*) buf_kop_only, (void*) NULL)) {
		if(BinTreeLookup(_canonicalParticipationCounts[topOrdinal][u][v], (foint)(void*) buf, (void*) &pcount))
		    ++*pcount;
		else {
		    pcount = Omalloc(sizeof(int));
		    *pcount = 1;
		    BinTreeInsert(_canonicalParticipationCounts[topOrdinal][u][v], (foint)(void*) buf, (foint)(void*) pcount);
		}
    #if 0
		int l,m;
		for(l=0;l<_k-1;l++) for(m=l+1;m<_k;m++) if(l!=i && m!=j && TinyGraphAreConnected(g,perm[l],perm[m])) {
		    int q=_orbitList[GintOrdinal][l];
		    int r=_orbitList[GintOrdinal][m];
		    int x=(int)perm[l];
		    int y=(int)perm[m];
		    if(x < y) { int tmp = x; x=y; y=tmp; }
		    if(r < q) { int tmp = q; q=r; r=tmp; }
		    // increase the count of octuplet (u,v, o,p, q,r, x,y)
		}
    #endif
	    }
        }
    }
}

// Given any canonical graphlet g, accumulate all submotifs of its canonical version. This is the
// fundamental pre-computation of the counts of (canonical node pair, canonical motif node pair)
// associations that's performed on the fly and then memoized for future use.
static void AccumulateCanonicalSubmotifs(int topOrdinal, TINY_GRAPH *g)
{
    static int depth;
    static Boolean initDone;
    static SET *seen; // size 2^B(k), *not* canonical but a specific set of nodes and edges in the *top* graphlet
    int i,j,l;
    if(!initDone) {
	assert(_Bk>0);
	seen = SetAlloc(_Bk);
	assert(_k>= 3 && _k <= 8);
	initDone = true;
    }

    int Gint = TinyGraph2Int(g,_k);
    if(depth==0){
	int GintOrdinal = _K[Gint];
	if(Gint != _canonList[GintOrdinal])
	    Fatal("AccumulateCanonicalSubmotifs can only initially be called with a canonical, but ord %d = %d != %d",
		GintOrdinal, _canonList[GintOrdinal], Gint);
	assert(GintOrdinal == topOrdinal);
	SetReset(seen);
    }
    if(SetIn(seen,Gint)) return;
    SetAdd(seen,Gint);
    SubmotifIncrementCanonicalPairCounts(topOrdinal, Gint, g);

    // Now go about deleting edges recursively.
    for(i=1; i<_k; i++)for(j=0;j<i;j++)
    {
	if(TinyGraphAreConnected(g,i,j)) // if it's an edge, delete it.
	{
	    TinyGraphDisconnect(g,i,j);
	    if(TinyGraphDFSConnected(g,0)) {
		++depth;
		AccumulateCanonicalSubmotifs(topOrdinal, g);
		--depth;
	    }
	    TinyGraphConnect(g,i,j);
	}
    }
}

/* This is called from ProcessGraphlet: a whole Varray of nodes from a sampled graphlet. Our job here is to
** accumulate the association counts for each pair of nodes, using the memoized counts from canonical graphlets
** computed above.
*/
void AccumulateGraphletParticipationCounts(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, int Gint, int GintOrdinal)
{
    static int ramCheck;
    if(++ramCheck % 100000 == 0) CheckRAMusage();
    int i,j;
    char perm[MAX_K];
    if(_canonicalParticipationCounts[GintOrdinal][1][0] == NULL) { // [1][0] since [0][0] will always be NULL
	for(i=1;i<_k;i++) for(j=0;j<i;j++)
	    _canonicalParticipationCounts[GintOrdinal][i][j] =
		BinTreeAlloc((pCmpFcn) strcmp, (pFointCopyFcn) strdup, (pFointFreeFcn) free, NULL, NULL);
	static TINY_GRAPH *canonical;
	if (!canonical) canonical = TinyGraphAlloc(_k);
	Int2TinyGraph(canonical, _canonList[GintOrdinal]);
	AccumulateCanonicalSubmotifs(GintOrdinal, canonical);
    }
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);

    double totalWeight = 0;
    for(i=0;i<_k;i++) totalWeight += 1.0/G->degree[Varray[i]];

    for(i=1;i<_k;i++) for(j=0;j<i;j++) // loop through all pairs of nodes in the *canonical* version of the graphlet
    {
	int g_u=perm[i], g_v=perm[j]; // u and v in the non-canonical motif g that is induced from G
	int G_u=Varray[g_u], G_v=Varray[g_v]; // u and v in the BIG input graph G.
#if PARANOID_ASSERTS
        assert(!TinyGraphAreConnected(g,g_u,g_v) == !GraphAreConnected(G,G_u,G_v));
#endif
	// Unlike the comment (+++) above, here we need info an *all* pairs of nodes in G that belong to this graphlet,
	// so we do not filter on the node pair being unconnected.
	if(G_u < G_v) { int tmp = G_u; G_u=G_v; G_v=tmp; } // for consistency in accessing _PredictGraph
	if(g_u < g_v) { int tmp = g_u; g_u=g_v; g_v=tmp; } // lower triangle of g
	_TraverseCanonicalPairs_G_u = G_u; _TraverseCanonicalPairs_G_v = G_v;
	_TraverseCanonicalPairs_Varray = Varray; _TraverseCanonicalPairs_perm = perm;
#if COUNT_uv_only && DEG_WEIGHTED_COUNTS
	_TraverseCanonicalPairs_deg_weight = totalWeight - 1.0/G->degree[G_u] - 1.0/G->degree[G_v];
#else
	_TraverseCanonicalPairs_deg_weight = 1;
#endif
	BinTreeTraverse(_canonicalParticipationCounts[GintOrdinal][i][j], TraverseCanonicalPairs);
#if 0
	int l,m;
	for(l=1;l<_k;l++) for(m=0;m<l;m++) if(l!=i && m!=j && TinyGraphAreConnected(g,perm[l],perm[m])) {
	    int q=perm[l], r=perm[m];
	    int x=Varray[q], y=Varray[r];
	    assert(GraphAreConnected(G,x,y));
	    if(x < y) { int tmp = x; x=y; y=tmp; }
	    if(q < r) { int tmp = q; q=r; r=tmp; }
	    PrintNodePairSorted(G_u,':',G_v);
	    printf(" %d %d:%d %d:%d", GraphAreConnected(G,G_u,G_v),g_u,g_v,q,r);
	    PrintNodePairSorted(x,':',y);
	    putchar('\n');
	}
#endif
    }

#if 0 //PREDICT_USE_TREE && have _THREADS
    Boolean debug = false;
    if(_child && _flushCounts) {
	fprintf(stderr, "Flushing child %d\n", getpid());
	_flushCounts = false;
	int G_u, G_v;
	for(G_u=1;G_u<G->n;G_u++) for(G_v=0;G_v<G_u;G_v++) {
	    assert(G_u > G_v);
	    if(PREDICT_GRAPH_NON_EMPTY(G_u,G_v)) {
		printf("%s %d", PrintNodePairSorted(G_u,':',G_v), GraphAreConnected(G,G_u,G_v));
		if(debug) {
		    _supportNodeNames=true;
		    fprintf(stderr, "CHILD has %d entries for %s %d", _PredictGraph[G_u][G_v]->n,
			PrintNodePairSorted(G_u,':',G_v), GraphAreConnected(G,G_u,G_v));
		    _supportNodeNames=false;
		}
		BinTreeTraverse(_PredictGraph[G_u][G_v], TraverseNodePairCounts);
		if(debug) fprintf(stderr, "\n");
		puts("");
		BinTreeFree(_PredictGraph[G_u][G_v]);
		_PredictGraph[G_u][G_v]=NULL;
	    }
	}
    }
#endif
}


int PredictMerge(GRAPH *G)
{
    assert(_outputMode == predict_merge);
    Predict_Init(G);
    if(isatty(fileno(stdin)))
	Warning("-mq (predict_merge) takes input only on stdin, which is currently a terminal. Press ^D or ^C to exit");
    assert(_JOBS==1); // we only read from standard input, so threads make no sense.
    char line[MAX_ORBITS * BUFSIZ];
    int lineNum = 0;
    while(fgets(line, sizeof(line), stdin) && !_memUsageAlarm) {
	Predict_ProcessLine(G, line);
	++lineNum;
    }
    assert(!ferror(stdin));
    assert(feof(stdin));
    PredictFlushAllCounts(G);
    return 0;
}
