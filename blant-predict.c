#include "Oalloc.h" // "Once alloc", for allocating tons of little integers via pointers
#include "../blant.h"
#include "../blant-output.h"
#include "../blant-utils.h"
#include "blant-predict.h"
#include "bintree.h"
#include "rand48.h"
#include <math.h>
#include <signal.h>
#include <ctype.h>

/* General note: when associating any pair of nodes (u,v) in the large graph G, with "orbit pairs" in canonical
** graphlets, the number of possible combination is so incredibly enormous and sparse, that even hashing was
** too wasteful of memory. Thus I moved to storing the associations for *each* pair (u,v) in one binary tree.
** Thus for G with n nodes, there are potentially (n choose 2) distinct binary trees, each holding the list
** of canonical node pair associations and the respective count. Basically for a 10,000 node graph, and using
** k=8, there are 80,000 possible orbits--except I discovered orbit pairs aren't enough, you need to use the
** actual nodes. Thus for k=8 there are 12346 canonical graphlets each with 7*8 possible pairs, for a total
** of 12346*7*8 = 691,000 possible canonical node pairings... so the total number of possible *associations*
** between any (u,v) in G and any canonical node pair is n^2 * 691,000... and for n=10,000 nodes, that's
** 7 x 10^13 associations... each of which must have a *count* (potentially requiring 32 bits), for a total
** of 3e14 bytes. So not feasible, and even a hash table (which is sparsely stored to reduce collisions)
** was requiring >100GB to store these associations. A Binary tree is a bit slower but requires only about 5GB.
*/
#if PREDICT_USE_BINTREE
static BINTREE ***_PredictGraph; // lower-triangular matrix (ie., j<i, not i<j) of dictionary entries

// Allocate the NULL pointers for just the *rows* of the PredictGraph[i].
void Predict_Init(GRAPH *G) {
#if PREDICT_USE_BINTREE
    int i;
    // This just allocates the NULL pointers, no actual binary trees are created yet.
    _PredictGraph = Calloc(G->n, sizeof(BINTREE**)); // we won't use element 0 but still need n of them
    for(i=1; i<G->n; i++) _PredictGraph[i] = Calloc(i, sizeof(BINTREE*));
#endif
}


/* Called indirectly via BinTreeTraverse, which prints *all* the currently stored participation counts
** of a particular pair of nodes (u,v) in G. This function is called on each such participation count,
** and it simply prints the k:g:i:j canonical set, and its count, to stdout.
*/
Boolean TraverseNodePairCounts(foint key, foint data) {
    char *ID = key.v;
    int *pCount = data.v;
    printf("\t%s %d",ID, *pCount);
    return true;
}
#endif


#define PREDICT_GRAPH_NON_EMPTY(i,j) (_PredictGraph[i][j] && _PredictGraph[i][j]->root)
/* Loop across all pairs of nodes (u,v) in G, and for each pair, print on one line the pair of nodes,
** the edge truth, and all the participation counts.
*/
void PredictFlushAllCounts(GRAPH *G){
	if(_child) ualarm(0); // turn off _flushCounts alarm
	int i,j;
#if !PREDICT_USE_AWK
	for(i=1; i < G->n; i++) for(j=0; j<i; j++) {
	    if(PREDICT_GRAPH_NON_EMPTY(i,j))  // only output node pairs with non-zero counts
	    {
		printf("%s %d", PrintNodePairSorted(i,':',j), GraphAreConnected(G,i,j));
#if PREDICT_USE_BINTREE
		BinTreeTraverse(_PredictGraph[i][j], TraverseNodePairCounts);
#endif
		puts("");
	    }
	}
#endif
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
    if(line[strlen(line)-1] != '\n')
	Fatal("char line[%s] buffer not long enough while reading child line in -mp mode",sizeof(line));
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
	IncrementNodePairCount(MAX(G_u,G_v) , MIN(G_u,G_v), ID, count);
	*s1 = tmp;
	assert(*(s0=s1)=='\n' || *s0 == '\t');
    }
    assert(*s0 == '\n');
}

#if PREDICT_USE_BINTREE
/* Rather than explicitly enumerating the sub-motifs of every graphlet sampled from G, we do the recursive enumeration
** of all motifs under a graphlet only one---on the canonical graphlet. Then we *store* the association between the
** nodes of the canonical graphlet, and the motif associations. The associations get stored in a binary tree, one
** for each pair of canonical nodes in each canonical graphlet. Then, when we sample a graphlet from G, we simple
** determine what canonical it is using our standard BLANT lookup table _K, and use the perm[] matrix to transfer
** the canonical associations to the nodes of the graphlet sampled from G. This function is called in that process,
** and the pair of nodes (u,v) in G must be global since they're set below and used here during the travelsal.
*/
static int _PredictGraph_u, _PredictGraph_v; // indices into _PredictGraph[u][v] since we don't know them during traversal
static Boolean TraverseCanonicalPairs(foint key, foint data) {
    // ID is a *string* of the form k:g:i:j, where g is the Ordinal of canonical graphlet g, and i and j are nodes in g
    char *ID = key.v;
    int *pCanonicalCount = data.v;
    return IncrementNodePairCount(_PredictGraph_u, _PredictGraph_v, ID, *pCanonicalCount);
}

// Given a pair of nodes (u,v) in G and an association ID, increment the (u,v) count of that association by the count of the ID.
Boolean IncrementNodePairCount(int G_u, int G_v, char *ID, int canonicalCount)
{
    int *pUVcount;
    if(G_u<G_v) { int tmp=G_u; G_u=G_v;G_v=tmp;}
    if(_PredictGraph[G_u][G_v] == NULL) _PredictGraph[G_u][G_v] =
	BinTreeAlloc(unbalanced, (pCmpFcn) strcmp, (pFointCopyFcn) strdup, (pFointFreeFcn) free, NULL, NULL);
    if(BinTreeLookup(_PredictGraph[G_u][G_v], (foint)(void*) ID, (void*) &pUVcount))
	*pUVcount += canonicalCount;
    else {
	pUVcount = Omalloc(sizeof(int));
	*pUVcount = canonicalCount;
	BinTreeInsert(_PredictGraph[G_u][G_v], (foint)(void*) ID, (foint)(void*) pUVcount);
#if PARANOID_ASSERTS
	assert(BinTreeLookup(_PredictGraph[G_u][G_v], (foint)(void*) ID, (void*) &pUVcount)
	    && *pUVcount == canonicalCount);
#endif
    }
    //printf("\t %s %d(%d)",ID, canonicalCount, *pUVcount);
    return true;
#endif
}


// a slightly more compact internal representation of the char*ID above, because sometimes we want to see and
// manipulate the values of g, i, and j, and it's cumbersome to extract them from a string (especially g).
typedef struct {
    short g; // the canonical ordinal (fits into short int)
    char i,j; // the two canonical nodes (0 through k-1)
} MOTIF_NODE_PAIR; // the count will be the "info" member of the binary tree node.

// ensure nodes i and j are *actually* disconnected in the motif
Boolean AssertMotifPairDisconnected(MOTIF_NODE_PAIR *op) {
    int GintOrdinal = op->g;
    int Gint = _canonList[GintOrdinal];
    static TINY_GRAPH *g;
    if(!g) g = TinyGraphAlloc(_k);  // only allocate this TinyGraph once, and re-use it each time.
    Int2TinyGraph(g, Gint);
    assert(!TinyGraphAreConnected(g,op->i,op->j));
    return true;
}

/*
** For each canonical graphlet, and for each pair of nodes in the graphlet, we maintain a dictionary (as a binary tree)
** of all the participation counts for all motif node-pairs that that canonical graphlet node pair participates in,
** across all sub-motifs of that canonical graphlet.
** Note that WE DO NOT KNOW HOW TO PROPERTLY NORMALIZE YET.... so for now do the easiest thing, which is count them all.
*/
static BINTREE *_canonicalParticipationCounts[MAX_CANONICALS][MAX_K][MAX_K];

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
	int u=(int)perm[i], v=(int)perm[j]; // u and v in the *current* motif
	// the association is between a node *pair* in the canonical top graphlet, and a node *pair* of the
	// motif that they are participating in. Since the pair is undirected, we need to choose a unique order
	// to use an a key, so we simply sort the pair.
	// And since we want lower triangle, row > column, always, so i=1...n-1, j=0...i-1
	if(u < v) { int tmp = u; u=v; v=tmp; }

        // We are trying to determine the frequency that a pair of nodes in the topOrdinal g have an edge based
        // on their being located at a pair of nodes in a motif. The frequency only makes sense if the underlying
        // edge between them can sometimes exist and sometimes not; but if TinyGraphAreConnected(u,v)=true, the
	// underlying edge is guaranteed to exist. Thus, we only want to check this node pair if the motif does NOT
	// have the edge. (Comment label: (+++))
        if(!TinyGraphAreConnected(g,u,v))
	{
	    MOTIF_NODE_PAIR op ={GintOrdinal,i,j};
#if PARANOID_ASSERTS
	    AssertMotifPairDisconnected(&op);
#endif
	    int *pcount;
	    char buf[BUFSIZ];
	    sprintf(buf,"%d:%d:%d:%d", _k, GintOrdinal,i,j);
	    if(_canonicalParticipationCounts[topOrdinal][i][j] == NULL)
		Fatal("oPC[%d][%d][%d] is null", topOrdinal,i,j);
	    if(BinTreeLookup(_canonicalParticipationCounts[topOrdinal][i][j], (foint)(void*) buf, (void*) &pcount)){
		++*pcount;
	    }
	    else {
		pcount = Omalloc(sizeof(int));
		*pcount = 1;
		BinTreeInsert(_canonicalParticipationCounts[topOrdinal][i][j], (foint)(void*) buf, (foint)(void*) pcount);
	    }
	    //printf("(topOrd,i,j)=(%d,%d,%d)   (Gord,u,v)=(%d,%d,%d) (o,p)=(%d,%d)=%d\n", topOrdinal,i,j, GintOrdinal,u,v, o,p,*pcount);
#if 0
	    int l,m;
	    for(l=0;l<_k-1;l++) for(m=l+1;m<_k;m++) if(l!=i && m!=j && TinyGraphAreConnected(g,perm[l],perm[m])) {
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

// Given any canonical graphlet g, accumulate all submotifs of its canonical version. This is the
// fundamental pre-computation of the counts of (canonical node pair, canonical motif node pair)
// associations that's performed on the fly and then memoized for future use.
void AccumulateCanonicalSubmotifs(int topOrdinal, TINY_GRAPH *g)
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
    assert(_child && !_flushCounts);
    fprintf(stderr, "Alarm in process %d\n", getpid());
    ualarm(0);
    _flushCounts = true;
}

/* This is called from ProcessGraphlet: a whole Varray of nodes from a sampled graphlet. Our job here is to
** accumulate the association counts for each pair of nodes, using the memoized counts from canonical graphlets
** computed above.
*/
void AccumulateGraphletParticipationCounts(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, int Gint, int GintOrdinal)
{
    int i,j;
    char perm[MAX_K];
    if(_canonicalParticipationCounts[GintOrdinal][1][0] == NULL) { // [1][0] since [0][0] will always be NULL
	for(i=1;i<_k;i++) for(j=0;j<i;j++)
	    _canonicalParticipationCounts[GintOrdinal][i][j] =
		BinTreeAlloc(unbalanced, (pCmpFcn) strcmp, (pFointCopyFcn) strdup, (pFointFreeFcn) free, NULL, NULL);
	static TINY_GRAPH *canonical;
	if (!canonical) canonical = TinyGraphAlloc(_k);
	Int2TinyGraph(canonical, _canonList[GintOrdinal]);
	AccumulateCanonicalSubmotifs(GintOrdinal, canonical);
    }
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);

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
#if PREDICT_USE_BINTREE
	_PredictGraph_u = G_u; _PredictGraph_v = G_v;
	BinTreeTraverse(_canonicalParticipationCounts[GintOrdinal][g_u][g_v], TraverseCanonicalPairs);
#else
	int l,m;
	for(l=0;l<_k-1;l++) for(m=l+1;m<_k;m++) if(l!=i && m!=j && TinyGraphAreConnected(g,perm[l],perm[m])) {
	    int x=Varray[(int)perm[l]];
	    int y=Varray[(int)perm[m]];
	    assert(GraphAreConnected(G,x,y));
	    if(x < y) { int tmp = x; x=y; y=tmp; }
	    if(r < q) { int tmp = q; q=r; r=tmp; }
	    PrintNodePairSorted(G_u,':',G_v);
	    printf(" %d %d:%d %d:%d", GraphAreConnected(G,G_u,G_v),o,p,q,r);
	    PrintNodePairSorted(x,':',y);
	    putchar('\n');
	}
#endif
    }

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
	signal(SIGALRM, AlarmHandler);
	ualarm(100000);  // (int)(2*_MAX_THREADS*drand48())); // on average, have one child per second send data to the parent process
    }
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
    while(fgets(line, sizeof(line), stdin)) {
	Predict_ProcessLine(G, line);
	++lineNum;
    }
    assert(!ferror(stdin));
    assert(feof(stdin));
    PredictFlushAllCounts(G);
    return 0;
}
