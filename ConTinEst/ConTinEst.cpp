/*

Implementation of the ConTinEst Class.  

Author : Nan Du (dunan@gatech.edu)

*/

#include "ConTinEst.h"

ConTinEst::ConTinEst(Graph *G, unsigned num_samples, unsigned num_rankings)
{
	m_RNG.SetState(0, 0);

	m_num_samples = num_samples;

	m_num_rankings = num_rankings;

	m_G_inverse = G;

	unsigned N = m_G_inverse->N;

	m_TableList = new vector<pair<float, unsigned> >**[m_num_samples];
	m_keys = new float**[m_num_samples];

	for (unsigned i = 0; i < m_num_samples; ++i)
	{
		m_TableList[i] = new vector<pair<float, unsigned> >*[m_num_rankings];
		m_keys[i] = new float*[m_num_rankings];

		for (unsigned j = 0; j < m_num_rankings; ++j)
		{
			m_TableList[i][j] = new vector<pair<float, unsigned> >[N];
			m_keys[i][j] = new float[N];
		}
	}
}

ConTinEst::~ConTinEst()
{
	unsigned N = m_G_inverse->N;
	
	for (unsigned i = 0; i < m_num_samples; ++i)
	{
		for (unsigned j = 0; j < m_num_rankings; ++j)
		{
			for(unsigned k = 0; k < N; ++ k)
			{
				vector<pair<float, unsigned> >().swap(m_TableList[i][j][k]);
			}
			delete[] m_TableList[i][j];
			delete[] m_keys[i][j];
		}
		delete[] m_TableList[i];
		delete[] m_keys[i];
	}

	delete[] m_TableList;
	delete[] m_keys;

	m_TableList = NULL;
	m_keys = NULL;
}

void ConTinEst::LeastElementListsSet(float *d, const pair<float, unsigned> *key_node_pairs, vector<pair<float, unsigned> > *lists)
{
	unsigned N = m_G_inverse->N;

	pair<float, unsigned> bound(-1,N);

	for(unsigned i = 0; i < N; ++ i)
	{
		unsigned vi = key_node_pairs[i].second;

		set<pair<float, unsigned> > key_to_node;
		map<unsigned, float> node_to_key;

		key_to_node.insert(make_pair(0, vi));
		node_to_key.insert(make_pair(vi, 0));

		while(!key_to_node.empty())
		{

			set<pair<float, unsigned> >::iterator itlow = key_to_node.lower_bound(bound);
			map<unsigned, float>::iterator itlow_con = node_to_key.find(itlow->second);

			float dk = itlow->first;
			unsigned vk = itlow->second;
			key_to_node.erase(itlow);
			node_to_key.erase(itlow_con);

			lists[vk].push_back(make_pair(dk, vi));
			d[vk] = dk;

			for (set<unsigned>::iterator c = m_G_inverse->nodes[vk].children.begin(); c != m_G_inverse->nodes[vk].children.end(); ++ c) {

				unsigned vj = *c;


				itlow_con = node_to_key.find(*c);

				float tmp = dk + m_G_inverse->edge_weight[vk][vj];
				if(itlow_con != node_to_key.end())
				{
					if(tmp < itlow_con->second)
					{
						itlow = key_to_node.find(make_pair(itlow_con->second, itlow_con->first));

						key_to_node.erase(itlow);
						key_to_node.insert(make_pair(tmp,itlow_con->first));
						itlow_con->second = tmp;

					}
				}else if(tmp < d[vj])
				{

					key_to_node.insert(make_pair(tmp, vj));
					node_to_key.insert(make_pair(vj, tmp));

				}


			}
			

		}


	}
}

void ConTinEst::GetLeastElementLists()
{
	unsigned N = m_G_inverse->N;

	pair<float, unsigned> *key_node_pairs = new pair<float, unsigned>[N];

	float *d = new float[N]; 
	
	for (unsigned i = 0; i < m_num_samples; ++ i) {

		m_G_inverse->SampleEdgeWeightWbl();

		for (unsigned j = 0; j < m_num_rankings; ++ j) {
		
			// initialize keys and distance d
			for (unsigned k = 0; k < N; ++ k) {
				float key = m_RNG.GetExponential(1.0);
				m_keys[i][j][k] = key;
				key_node_pairs[k].first = key;
				key_node_pairs[k].second = k;
			}

			sort(key_node_pairs, key_node_pairs + N);
			fill_n(d, N, 1e10);

			LeastElementListsSet(d, key_node_pairs, m_TableList[i][j]);
		}

	}

	delete[] d;
	delete[] key_node_pairs;	
}

float ConTinEst::EstimateNeighborhood(const set<unsigned>& sources, float T)
{

	float size = 0.0, avg = 0.0;
	unsigned N = m_G_inverse->N;

	if(!sources.empty())
	{
		for (unsigned i = 0; i < m_num_samples; ++ i) {

			size = 0.0;

			for (unsigned j = 0; j < m_num_rankings; ++ j) {

				float minRank = 1e10;
				pair<float, unsigned> tmp(T, N);
				
				for (set<unsigned>::const_iterator s = sources.begin(); s != sources.end(); ++ s) {
				
					vector<pair<float, unsigned> >::iterator idx = lower_bound(m_TableList[i][j][*s].begin(), m_TableList[i][j][*s].end(), tmp, greater<pair<float, unsigned> >());
				
					if (m_keys[i][j][idx->second] < minRank) {
						minRank = m_keys[i][j][idx->second];
					}
				
				}

				size += minRank;

			}

			avg += ((m_num_rankings - 1) / size);

		}

		return avg / m_num_samples;
	}else
	{
		return 0;
	}

}


set<unsigned> ConTinEst::LZGreedy(double T, unsigned K)
{
	unsigned N = m_G_inverse->N;

	set<unsigned> sources;

	vector<pair<double, unsigned> > marginal_gain;
	make_heap(marginal_gain.begin(), marginal_gain.end());
	
	for (unsigned i = 0; i < N; ++i)
	{
		set<unsigned> tmp;
		tmp.insert(i);

		marginal_gain.push_back(make_pair(EstimateNeighborhood(tmp, T), i));
		push_heap(marginal_gain.begin(), marginal_gain.end());

	}

	pair<double, unsigned> &top_max = marginal_gain.front();

	double total = top_max.first;
	sources.insert(top_max.second);

	pop_heap(marginal_gain.begin(), marginal_gain.end());
	marginal_gain.pop_back();

	bool *valid = new bool[N];

	while(sources.size() < K)
	{
		fill_n(valid, N, false);

		while(true)
		{
			top_max = marginal_gain.front();

			if (valid[top_max.second])	
			{
				sources.insert(top_max.second);
				total += top_max.first;
				pop_heap(marginal_gain.begin(), marginal_gain.end());
				marginal_gain.pop_back();
				break;
			}

			set<unsigned> tmp = sources;
			tmp.insert(top_max.second);

			double gain = EstimateNeighborhood(tmp, T) - total;



			top_max.first = (gain > 0) ? gain : 0;
			valid[top_max.second] = true;
			make_heap(marginal_gain.begin(), marginal_gain.end());

		}
	}

	// for(set<unsigned>::const_iterator s = sources.begin(); s != sources.end(); ++ s)
	// {
	// 	cout << *s << " ";
	// }
	// cout << endl;


	delete[] valid;

	return sources;
}

vector<set<unsigned> > ConTinEst::Optimize(const vector<double>& setT, const vector<unsigned>& setK)
{
	vector<set<unsigned> > sources;
	for(vector<unsigned>::const_iterator k = setK.begin(); k != setK.end(); ++ k)
	{
		for(vector<double>::const_iterator t = setT.begin(); t != setT.end(); ++ t)
		{
			sources.push_back(LZGreedy(*t, *k));
		}
	}
	
	return sources;	
}