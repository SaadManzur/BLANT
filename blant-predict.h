#ifndef BLANT_PREDICT_H
#define BLANT_PREDICT_H

// Method used to store the sparse relationship between node pairs in G, and Canonical Graphlet Node Pairs.
// BINTREE is actually the only one implemented--hashmap required tens to hundreds of GB while BT needs only a handful
// and is plenty fast enough.
#define PREDICT_USE_BINTREE 1
#define PREDICT_USE_AWK !(PREDICT_USE_BINTREE) // use awk for the associations
#if (PREDICT_USE_BINTREE+PREDICT_USE_AWK) != 1
#error can only choose one method of accumulating node pair associations
#endif

void Predict_Init(GRAPH *G);
void PredictFlushAllCounts(GRAPH *G);
void Predict_ProcessLine(GRAPH *G, char line[]);

void AccumulateGraphletParticipationCounts(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, int Gint, int GintOrdinal);
int PredictMerge(GRAPH *G);

#endif
