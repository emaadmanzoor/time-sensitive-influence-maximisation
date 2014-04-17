/*

Definition of the Naive Simulation Class.  

Author : Nan Du (dunan@gatech.edu)

*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Common.h"
#include "Graph.h"

using namespace std;

struct Node2Time{
	unsigned nodeID;
	double time;
	
	bool operator < (const Node2Time& n2t) const{
		if (time < n2t.time) {
			return true;
		}else {
			return false;
		}
	}
};

class Simulation{

private:

	Graph* G;

	SimpleRNG RNG;

private:

	void GenerateCascade(vector<Node2Time>& cascade, set<unsigned>& initialSet, double TimeHorizon, map<unsigned, unsigned>& infectedBy);

public:
	Simulation(Graph *G);
	~Simulation(){};

	double RandomSimulation(double T, set<unsigned>& initialSet, unsigned C);


};

#endif