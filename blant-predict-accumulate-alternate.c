/*
** NOTE 1/22/2021: Up until now, I think we've been doing it wrong, as follows: when counting *actual* L3s between nodes
** that actually have an edge between them in the original graph, it shouldn't matter whether the pair of nodes are in the
** training set (in which case the edge is observed), or the test set (in which case the edge is absent in the training
** set--which DOESN'T CHANGE THE L3 COUNT OF THESE TWO NODES!!!). Thus, if you compare the *distribution* of actual L3
** counts of true edges in the training set, to that of node pairs that are in the training set, the distributions are
** *identical*, as they should be.
**    However, this fact is *not* true of the counts associated with orbit pair 4:11:11; node pairs with observed
** edges have about 1.7x higher 4:11:11 count than node pairs whose edge has been moved to the test set. I think this
** is because I don't delete edges in the sampled graphlet *until* we start looking at its motifs. That's incorrect,
** and your long-ago original Kovacs algorithm was closer: you need to pick a pair of nodes in the sampled *graphlet*,
** and *then* sample the motifs... and if the pair of nodes is an edge, you should delete that edge BEFORE YOU START
** looking at the motifs underneath. Thus, another level of edge deletion has to occur (which only adds to pre-computation):
** you need to loop through all all pairs of nodes at the top-level graphlet, and delete the edge if it's their; and then
** recurse down the motifs but counting *ONLY* the node pair at the top. So in the recursion we need to keep track
** not only of the TopOrdinal, but the Top_u and Top_v who's orbit-pairs we're accumulating, and update the orbit-pairs
** of *only* that node pair as we recurse down. But then there needs to be (k choose 2) top-level calls to the recursive
** motif algorithm.
**
** NOTE 2: The above was 3pm. I wrote all of the below.... but it seems to make results *worse*. Go figure???
*/

// Given a topOrdinal and two of its nodes top_u, top_v, accumulate the orbit pair participation counts
// of *just* that pair of nodes (top_u,top_v).
void SubmotifIncrementCanonicalPairCounts(int topOrdinal, int top_u, int top_v, TINY_GRAPH *g)
{
    int Gint = TinyGraph2Int(g,_k);
    int GintOrdinal=_K[Gint];
    char perm[MAX_K], inversePerm[MAX_K];
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);
    InvertPerm(inversePerm, perm);

    // We are trying to determine the frequency that a pair of nodes in the topOrdinal have an edge based on their
    // being located at a pair of canonical nodes in a sub-motif. The frequency only makes sense if the underlying
    // edge between them can sometimes exist and sometimes not; but if TinyGraphAreConnected(u,v)=true, the edge
    // already exists in this motif (and thus also in the topOrdinal) and there's nothing to predict. Thus, we only
    // want to check this node pair if the motif does NOT have the edge. (Comment label: (+++))

    int i=inversePerm[top_u], j=inversePerm[top_v]; // i and j are now the canonical node IDs

    assert(!TinyGraphAreConnected(g,top_u,top_v));
    int o=_orbitList[GintOrdinal][i], p=_orbitList[GintOrdinal][j];
    // the association is between a *node pair* in the canonical top graphlet, and an *orbit pair* of the
    // motif that they are participating in. Since the pair is undirected, we need to choose a unique order
    // to use an a key, so we simply sort the pair.
    // And since we want lower triangle, row > column, always, so o=1...n-1, p=0...o-1
    if(top_u < top_v) { int tmp = top_u; top_u=top_v; top_v=tmp; }
    if(o < p) { int tmp = o; o=p; p=tmp; }
    assert(_canonicalParticipationCounts[topOrdinal][top_u][top_v]);
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
	if(BinTreeLookup(_canonicalParticipationCounts[topOrdinal][top_u][top_v], (foint)(void*) buf, (void*) &pcount))
	    ++*pcount;
	else {
	    pcount = Omalloc(sizeof(int));
	    *pcount = 1;
	    BinTreeInsert(_canonicalParticipationCounts[topOrdinal][top_u][top_v], (foint)(void*) buf, (foint)(void*) pcount);
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
	    // increase the count of octuplet (top_u,top_v, o,p, q,r, x,y)
	}
#endif
    }
}


static void Helper_AccumulateCanonicalSubmotifs(int topOrdinal, int top_u, int top_v, TINY_GRAPH *g)
{
    int i,j;
#if PARANOID_ASSERTS
    assert(!TinyGraphAreConnected(g,top_u,top_v));
#endif

    SubmotifIncrementCanonicalPairCounts(topOrdinal, top_u, top_v, g);

    // Now go about deleting edges recursively to enumerate all connected motifs underneath
    for(i=1; i<_k; i++)for(j=0;j<i;j++)
    {
	if(TinyGraphAreConnected(g,i,j)) // if it's an edge, delete it.
	{
	    TinyGraphDisconnect(g,i,j);
	    if(TinyGraphDFSConnected(g,0))
		Helper_AccumulateCanonicalSubmotifs(topOrdinal, top_u, top_v, g);
	    TinyGraphConnect(g,i,j);
	}
    }
}

// Given any canonical graphlet g, accumulate all submotifs of its canonical version. This is the
// fundamental pre-computation of the counts of (canonical node pair, canonical motif node pair)
// associations that's performed on the fly and then memoized for future use.

static void AccumulateCanonicalSubmotifs(int topOrdinal, TINY_GRAPH *g)
{
    int Gint = TinyGraph2Int(g,_k);
    int GintOrdinal = _K[Gint];
    if(Gint != _canonList[GintOrdinal])
	Fatal("AccumulateCanonicalSubmotifs can only initially be called with a canonical, but ord %d = %d != %d",
	    GintOrdinal, _canonList[GintOrdinal], Gint);
    assert(GintOrdinal == topOrdinal);

    int i,j;
    for(i=1; i<_k; i++)for(j=0;j<i;j++)
    {
	// As described above, we need to delete the edge BEFORE starting to count the submotifs.
	Boolean edge = TinyGraphAreConnected(g,i,j);
	if(!edge) Helper_AccumulateCanonicalSubmotifs(topOrdinal, i, j, g);
	else {
	    TinyGraphDisconnect(g,i,j);
	    if(TinyGraphDFSConnected(g,0)) // if removing the edge didn't disconnect the graphlet...
		Helper_AccumulateCanonicalSubmotifs(topOrdinal, i, j, g); // exactly the same call as above
	    TinyGraphConnect(g,i,j); // reconnect the edge
	}
    }
}

