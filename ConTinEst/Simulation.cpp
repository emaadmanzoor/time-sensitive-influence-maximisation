/*

Implementation of the Naive Simulation Class.  

Author : Nan Du (dunan@gatech.edu)

*/

#include "Simulation.h"

Simulation::Simulation(Graph *G_wbl)
{
	G = G_wbl;

	RNG.SetState(0, 0);
}

void Simulation::GenerateCascade(vector<Node2Time>& cascade, set<unsigned>& initialSet, double TimeHorizon, map<unsigned, unsigned>& infectedBy)
{
	double GlobalTime = 0.0;
	
	vector<Node2Time> visitedNodes;
	
	for (set<unsigned>::iterator u = initialSet.begin(); u != initialSet.end(); ++ u) {
		Node2Time n2t;
		n2t.nodeID = *u;
		n2t.time = GlobalTime;
		visitedNodes.push_back(n2t);
	}
	
	while (true) {
		
		sort(visitedNodes.begin(), visitedNodes.end());
		
		
		map<unsigned, unsigned> visitedNodes2Idx;
		for (unsigned i = 0; i < visitedNodes.size(); ++ i) {
			visitedNodes2Idx.insert(make_pair(visitedNodes[i].nodeID, i));
		}
		
		unsigned currentID = visitedNodes[0].nodeID;
		GlobalTime = visitedNodes[0].time;
		if(GlobalTime >= TimeHorizon)
		{
			break;
		}
		
		cascade.push_back(visitedNodes[0]);
		
		for (set<unsigned>::iterator u = G->nodes[currentID].children.begin(); u != G->nodes[currentID].children.end(); ++ u) {
			
			
			if((infectedBy.find(currentID) != infectedBy.end()) && (infectedBy[currentID] == *u))
			{
				continue;
			}
			
			double deltaT = RNG.GetWeibull(G->edge_parameter[currentID][*u].shape,G->edge_parameter[currentID][*u].scale);
					
			double t1 = GlobalTime + deltaT;
			
			map<unsigned, unsigned>::iterator m = visitedNodes2Idx.find(*u);
			
			if(m != visitedNodes2Idx.end())
			{
				double t2 = visitedNodes[m->second].time;
				if((t2 != TimeHorizon) && (t2 > t1))
				{
					visitedNodes[m->second].time = t1;
					infectedBy[*u] = currentID;
				}
			}else {
				Node2Time n2t;
				n2t.nodeID = *u;
				n2t.time = t1;
				
				visitedNodes2Idx.insert(make_pair(*u, visitedNodes.size()));
				visitedNodes.push_back(n2t);
				infectedBy.insert(make_pair(*u, currentID));
			}
			
		}
		
		visitedNodes[0].time = TimeHorizon;
		
	}

}

double Simulation::RandomSimulation(double T, set<unsigned>& initialSet, unsigned C)
{
	vector<Node2Time> cascade;
	
	vector<double> Node2Infected(G->N, 0);
	
	map<unsigned, unsigned> infectedBy;
	
	for (unsigned c = 0; c < C; ++ c) {
		
		cascade.clear();
		
		infectedBy.clear();
		
		if (initialSet.size() > 0) {
			
			GenerateCascade(cascade, initialSet, T, infectedBy);
			
			for (int i = 0; i < cascade.size(); ++ i) {
				if (cascade[i].time <= T) {
					Node2Infected[cascade[i].nodeID] += 1;
				}
			}
		}
		
		
	}
	
	for (unsigned i = 0; i < Node2Infected.size(); ++ i) {
		Node2Infected[i] = Node2Infected[i] / C;
	}

	double sum = 0;
	for (vector<double>::const_iterator j = Node2Infected.begin(); j != Node2Infected.end(); ++j)
	{
		sum += *j;
	}

	return sum;
	
}