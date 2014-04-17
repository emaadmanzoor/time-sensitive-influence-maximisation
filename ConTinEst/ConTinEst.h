/*

Definition of the ConTinEst Class.  

Author : Nan Du (dunan@gatech.edu)

*/

#ifndef CONTINEST_H
#define CONTINEST_H

#include "Common.h"
#include "Graph.h"

using namespace std;

class ConTinEst{

private:

	SimpleRNG m_RNG;

	unsigned m_num_samples;

	unsigned m_num_rankings;

	Graph *m_G_inverse;

	vector<pair<float, unsigned> >*** m_TableList;

	float *** m_keys;


private:

	void LeastElementListsSet(float *d, const pair<float, unsigned> *key_node_pairs, vector<pair<float, unsigned> > *lists);

	set<unsigned> LZGreedy(double T, unsigned K);

public:

	ConTinEst(Graph *G, unsigned num_samples, unsigned num_rankings);
	~ConTinEst();

	void GetLeastElementLists();

	float EstimateNeighborhood(const set<unsigned>& sources, float T);

	vector<set<unsigned> > Optimize(const vector<double>& setT, const vector<unsigned>& setK);

};


#endif