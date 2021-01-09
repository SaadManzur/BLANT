#ifndef BLANT_PREDICT_H
#define BLANT_PREDICT_H

void Predict_Init(GRAPH *G);
void PredictFlushAllCounts(GRAPH *G);
void Predict_ProcessLine(GRAPH *G, char line[]);

void AccumulateGraphletParticipationCounts(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, int Gint, int GintOrdinal);
int PredictMerge(GRAPH *G);

#endif
