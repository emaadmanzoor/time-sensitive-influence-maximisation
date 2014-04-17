/*

Utility Class for Unit Test.  

Author : Nan Du (dunan@gatech.edu)

*/

#include "TestModules.h"

 void TestModules::TestConTinEst(string filename, unsigned int N)
 {
 	cout << filename << endl;

	Graph G(filename, N);

	G.LoadWeibullFormatNetwork(",", true);

	unsigned num_samples = 10000, num_labels = 5;

	ConTinEst continest(&G, num_samples, num_labels);

	cout <<"Get all least-label lists : " << num_samples << " sets of transmission times; " << num_labels << " sets of random labels; ";
	continest.GetLeastElementLists();
	cout <<"done" << endl << endl;


	Graph G1(filename, N);
	G1.LoadWeibullFormatNetwork(",", false);

	pair<unsigned, unsigned> result = G1.MaximumOutDegree();

	unsigned nodeID = result.first;
	unsigned degree = result.second;

	cout << "node " << nodeID << " has the largest out-degree " << degree << endl;

	set<unsigned> sources;
	sources.insert(nodeID);

	unsigned C = 10000;

	for(unsigned T = 1; T <= 10; ++ T)
	{
		cout << "Estimate Influence T = " << T;
		double estimated_influence = continest.EstimateNeighborhood(sources, T);
		cout << " done" << endl << endl;

		
		cout << "Simulation Check T = " << T << " " << " with " << C << " samples" << endl;
		Simulation MC(&G1);

		cout << "ConTinEst : " << estimated_influence << endl << "Simulation : " << MC.RandomSimulation(T, sources, C) << endl << endl;
	}

	vector<double> set_T;
	vector<unsigned> set_K;

	set_T.push_back(10);
	set_K.push_back(10);

	cout <<"Influence Maximization by selecting 10 nodes with T = 10 ";

	vector<set<unsigned> > tables = continest.Optimize(set_T, set_K);

	cout << "done" << endl;

	cout << "selected sources : " ;

	for(vector<set<unsigned> >::const_iterator m = tables.begin(); m != tables.end(); ++ m)
	{
		for(set<unsigned>::const_iterator u = m->begin(); u != m->end(); ++ u)
		{
			cout << *u << ";";
		}
		cout << endl;
	}


 }