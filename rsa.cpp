#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <algorithm>
#include "rsa.h"
// ----- KSP ---------------------
#include <limits>
#include <set>
#include <map>
#include <queue>
#include <vector>
#include "GraphElements.h"
#include "Graph.h"
#include "DijkstraShortestPathAlg.h"
#include "YenTopKShortestPathsAlg.h"
// -------------------------------
// ----- QueryPerformanceCounter -----
#include <windows.h>
#include <Mmsystem.h>
#pragma comment(lib, "Winmm.lib")
// -----------------------------------

// ----- function definition -----
void ShowTwoDimArray(int, int, int**);  // for test
void OutputLinkTopoFile(Topology , string);  // for Routing::KShortestPath(int, Topology).
int** CreateTwoDimArray(int row, int col, int value);  // for create 2-dimension int array.
EdgeTransformation SetET(Topology);  // to initialize the EdgeTransformation struct.
void ShowVector(string, vector<int>);
// -------------------------------

// ----- function for class Routing ------------
Routing::Routing()  // constructor
{
    // nothing
}

Routing::~Routing()  // destructor
{
    // delete path[N][N][K] and piarK[N][N]
	int N = Node_num;
    try
    {
		for(int i=0; i<N; i++)
    	{
			//for(int j=0; j<N; j++)
	    		//delete [] path[i][j];
			delete [] path[i];
			delete [] pairK[i];
			//delete [] path_set[i];
			//delete [] tmp[i];
		}
		delete [] path;
		delete [] pairK;
		//delete [] path_set;
		//delete [] tmp;
	}
	catch(const char* message)
	{
		cout<<"Destruct class Routing failed!"<<endl;
		exit(0);
	}
}

int Routing::Get_K_value()
{
	return this->Kpath;
}

int Routing::Get_Path_Hop(PhysicalPath p)
{
	return p.node.size()-1;
}

/*
void Routing::AllShortestPath(int K, Topology G)
{
	// set variable
	this->Kpath = K;
	this->Node_num = G.node_num;
	// output link-based topology file
	OutputLinkTopoFile(G, "topology.txt");
	// initialize the Path data structure: path[N][N][K]
	int N = G.node_num;
	path_set = new vector<PhysicalPath>*[N];
	path_time = new double*[N];
	for(int i=0; i<N; i++)
	{
		path_set[i] = new vector<PhysicalPath>[N];
		path_time[i] = new double[N];
	}
	tmp = new vector<double> *[N];
	for(int i=0; i<N; i++)
	{
		tmp[i] = new vector<double>[N];
	}
	// couting time (QueryPerformanceCounter)
	LARGE_INTEGER  large_interger;
	QueryPerformanceFrequency(&large_interger);
	double dff = large_interger.QuadPart;
	__int64 tp1;
	// running KSP
	Graph my_graph("input/topology.txt");
	for(int i=0; i<N; i++)  // for each (s,d) pair
	{
		for(int j=0; j<N; j++)
		{
			//if(i==j)
				//continue;

			// start to count time
			QueryPerformanceCounter(&large_interger);
			tp1 = large_interger.QuadPart;
			//double tp = clock();

			// calculate each SP in K value
			YenTopKShortestPathsAlg yenAlg(my_graph, my_graph.get_vertex(i), my_graph.get_vertex(j));
			while(yenAlg.has_next())
			{
				QueryPerformanceCounter(&large_interger);//
				__int64 tp2 = large_interger.QuadPart;//
				BasePath* p = yenAlg.next();
				PhysicalPath tmp_p;
				//p->PrintOut(cout);  // print out the path
				tmp_p.len = p->Weight();
				tmp_p.node.reserve(p->length());				
				for(std::vector<BaseVertex*>::const_iterator pos=p->m_vtVertexList.begin(); pos!=p->m_vtVertexList.end();++pos)
				{
					tmp_p.node.push_back((*pos)->getID());
				}
				path_set[i][j].push_back(tmp_p);
				//ShowPath(tmp_p);
				QueryPerformanceCounter(&large_interger);//
				tp2 = large_interger.QuadPart - tp2;//
				//cout<<current_k++<<": "<<tp2*1000/dff<<endl;//
				if(tmp[i][j].empty())
					tmp[i][j].push_back(tp2*1000/dff);
				else
				{
					double pre = tmp[i][j].back();
					tmp[i][j].push_back(pre+tp2*1000/dff);
				}
			}
			//path_time[i][j] = clock()-tp;
			// end to count time
			QueryPerformanceCounter(&large_interger);
			tp1 = large_interger.QuadPart - tp1;
			path_time[i][j] = tp1*1000/dff;
			//cout<<path_time[i][j]<<endl;
			//exit(0);
		}
	}
	//for(int i=0; i<N; i++)
		//for(int j=0; j<N; j++)
			//printf("[%2d][%2d]:\t%4f\tms\n",i,j,path_time[i][j]);
	//exit(0);
}
*/

void Routing::KShortestPath(int K, Topology G)
{
	// set variable
	this->Kpath = K;
	this->Node_num = G.node_num;
	// output link-based topology file
	OutputLinkTopoFile(G, "topology.txt");
	// initialize the Path data structure: path[N][N][K]
	int N = G.node_num;
	path = new vector<PhysicalPath>*[N];
	for(int i=0; i<N; i++)
	{
		path[i] =  new vector<PhysicalPath>[N];
	}
	// initialize the record of each path's K value: pairK[N][N]
	pairK = new int*[N];
	for(int i=0; i<N; i++)
	{
		pairK[i] = new int[N];
	}
	// running KSP
	Graph my_graph("input/topology.txt");
	for(int i=0; i<N; i++)  // for each (s,d) pair
	{
		for(int j=0; j<N; j++)
		{
			// calculate each SP in K value
			YenTopKShortestPathsAlg yenAlg(my_graph, my_graph.get_vertex(i), my_graph.get_vertex(j));
			int current_k = 0;
			while(current_k<K && yenAlg.has_next())
			{
				BasePath* p = yenAlg.next();
				//p->PrintOut(cout);  // print out the path
				PhysicalPath get_p;
				get_p.len = p->Weight();
				get_p.node.reserve(p->length());			
				for(std::vector<BaseVertex*>::const_iterator pos=p->m_vtVertexList.begin(); pos!=p->m_vtVertexList.end();++pos)
				{
					get_p.node.push_back((*pos)->getID());
				}
				this->path[i][j].push_back(get_p);
				++current_k;
			}
			pairK[i][j] = current_k;
		}
	}	
}

PhysicalPath* Routing::ReturnKPath(int S, int D)
{
	return &path[S][D][0];
	//return &path_set[S][D][0];
}

PhysicalPath Routing::ReturnSinglePath(int S, int D, int k)
{
	return this->path[S][D][k];
	//return this->path_set[S][D][k];
}

int Routing::GetPathNum(int S, int D)
{
	//return path_set[S][D].size();
	return path[S][D].size();
}

/*
void Routing::ShowSingleKPath(int sour, int dest, int k)
{
	if(k >= Kpath || k >= pairK[sour][dest])
		return;
	// show the path[sour][dest][k]
	cout<<path[sour][dest][k].node[0];
	for(int i=1; i<path[sour][dest][k].node.size(); i++)
		cout<<"-"<<path[sour][dest][k].node[i];
	cout<<endl;
}
*/

void Routing::ShowPath(PhysicalPath p)
{
	cout<<p.node[0];
	for(int i=1; i<p.node.size(); i++)
		cout<<"-"<<p.node[i];
	cout<<endl;
}
// ---------------------------------------------

// ----- function for class cbrsa --------------
cbrsa::cbrsa(Topology mainG, int Fmax, int ReqMaxSize)
{
	this->Fmax = Fmax;
	this->ReqMaxSize = ReqMaxSize;
	// initialize pool
	pool.num = 1;
	// create pool.block[Fmax][2]
	pool.block = CreateTwoDimArray(Fmax, 2, 0);
    pool.block[0][0] = 0;
    pool.block[0][1] = Fmax-1;
	pool.sum = Fmax;
	// assign Topology G
	G.node_num = mainG.node_num;
	G.graph.assign(mainG.graph.begin(), mainG.graph.end());
	// initialize et
    this->et = SetET(mainG);
    	//ShowTwoDimArray(et.node_num, et.node_num, et.eid);
	// initalize working and backup cluster
    Cw = new Cluster[Fmax];
	for(int i=0; i<Fmax; i++)
	{
		Cw[i].weight = new int[et.edge_num];
		Cw[i].SizeBucket.reserve(ReqMaxSize+1);
		Cw[i].SizeBucket.assign(ReqMaxSize+1,0);
	}
	sidx_Cw_seq.clear();  // empty;
	sidx_Cw_seq.reserve(Fmax);
	mf_Cw_seq = new vector< pair<int,int> >[ReqMaxSize+1];
	Cw_size = new vector<int>[ReqMaxSize+1];
	for(int i=0; i<ReqMaxSize+1; i++)
	{
		mf_Cw_seq[i].reserve(Fmax);
		Cw_size[i].reserve(Fmax);
	}
	waiting_queue.clear();
	req_info.clear();
	wPath_set.clear();
	for(int i=Fmax-1; i>=0; i--)  // push the candidate working cluster id into candidate_wcluster_id.
		candidate_wcluster_id.push_back(i);

	// ------ initial analysis date ----------
	r_count = 0;
	size_block = 0;
	total_Cw_num = 0;
	sum_Cw_size = 0;
	sum_ctn_cw = 0;
	total_super_num = 0;
}

cbrsa::~cbrsa()
{
	// delete pool[Fmax][2], Cw[Fmax], Cw[F].weight[et.edge_num], et.eid[N][N], mf_Cw_seq[Max_Req_Size]
	int F = this->Fmax;
	int RmaxSize = this->ReqMaxSize;
	int N = et.node_num;
	for(int i=0; i<F; i++)
	{
		delete [] pool.block[i];
		delete [] Cw[i].weight;
	}
	for(int i=0; i<N; i++)
		delete [] et.eid[i];
	delete [] pool.block;
	delete [] et.eid;
	delete [] mf_Cw_seq;
	delete [] Cw_size;
	delete [] Cw;
}

// ----- get/set value -----
int cbrsa::Get_Cw_sidx(int cid)
{
	return Cw[cid].sidx;
}

PhysicalPath cbrsa::GetwPath(int req_id)
{
	return req_info.find(req_id)->second.wp;
}

void cbrsa::RemoveCwFromSidxCwSeq(int cid)
{
	sidx_Cw_seq.erase(std::remove(sidx_Cw_seq.begin(), sidx_Cw_seq.end(), cid), sidx_Cw_seq.end());
}

void cbrsa::AddCwIntoMfCwSeq(int cid, int size)
{
	for(int i=1; i<=size; i++)
	{
		mf_Cw_seq[i].push_back(std::pair<int,int>(size-i,cid));
		std::sort(mf_Cw_seq[i].begin(),mf_Cw_seq[i].end());
	}
}

void cbrsa::RemoveCwFromMfCwSeq(int cid, int size)
{
	for(int i=1; i<=size; i++)
		mf_Cw_seq[i].erase(std::remove(mf_Cw_seq[i].begin(), mf_Cw_seq[i].end(), pair<int,int>(size-i,cid)), mf_Cw_seq[i].end());
}

void cbrsa::RemoveCwFromCwSize(int cid, int size)
{
	Cw_size[size].erase(std::remove(Cw_size[size].begin(), Cw_size[size].end(), cid), Cw_size[size].end());
}

void cbrsa::AddIntoWaitingQueue(int req_id, int size, int cid, PhysicalPath p)
{
	// create request
	RequestStruct req;
	req.size = size;
	req.MultiCw = false;
	req.wcid_set.push_back(cid);
	req.size_set.push_back(size);
	req.wp.len = p.len;
	req.wp.node.assign(p.node.begin(), p.node.end());
	// add into waiting_queue[]
	waiting_queue.insert(std::make_pair(req_id,req));
}

void cbrsa::AddIntoWaitingQueue(int req_id, int size, vector<int> mcid, vector<int> msize, PhysicalPath p)
{
	//cout<<"req "<<req_id<<" use super-cluser: "<<endl;
	//ShowVector("super-cluster[]", mcid);
	//ShowVector("size_set[]", msize);

	// create request
	RequestStruct req;
	req.size = size;
	req.MultiCw = false;
	req.wcid_set.assign(mcid.begin(), mcid.end());
	req.size_set.assign(msize.begin(), msize.end());
	req.wp.len = p.len;
	req.wp.node.assign(p.node.begin(), p.node.end());
	// add into waiting_queue[]
	waiting_queue.insert(std::make_pair(req_id,req));
}

// ---------------------

// ----- method -----
int cbrsa::FindCwToAssign_KSP_SP(int req_id, int size, int Kpath, PhysicalPath p[])
{
	if(!sidx_Cw_seq.empty())
	{
		for(int k=0; k<Kpath; k++)  // for each candidate path
		{
			// find the existing cluster to assign req
			for(std::vector<int>::size_type i=0; i!=sidx_Cw_seq.size(); i++)
			{
				int cluster_id = sidx_Cw_seq[i];
				if(Cw[cluster_id].size < size)  // the cluster size < req's size
					continue;
				if(CheckCwAssign(cluster_id, size, p[k]))
				{
					// p[k] is the selected path for req., and add p[k] into Path_set.
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[k]));
					return cluster_id;
				}
			}
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
	}
	return ci;
}

int cbrsa::FindCwToAssign_KSP_FF(int req_id, int size, int Kpath, PhysicalPath p[])
{
	if(!sidx_Cw_seq.empty())
	{
		for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
		{
			int cluster_id = sidx_Cw_seq[i];
			for(int k=0; k<Kpath; k++)  // for each candidate path
			{
				if(Cw[cluster_id].size < size)  // the cluster size < req's size
					continue;
				if(CheckCwAssign(cluster_id, size, p[k]))
				{
					// assign this path p[k] to req, and jump out the k-path for-loop
					// p[k] is the selected path for req., and add p[k] into Path_set.
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[k]));
					return cluster_id;
				}
			}
		}
	}

	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
	}
	return ci;
}

int cbrsa::FindCwToAssign_KSP_MF(int req_id, int size, int Kpath, PhysicalPath p[])
{
	if(!mf_Cw_seq[size].empty())
	{
		// get cluster id
		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
		{			
			int cluster_id = it->second;
			for(int k=0; k<Kpath; k++)  // for each candidate path
			{
				if(Cw[cluster_id].size < size)  // the cluster size < req's size
					continue;
				if(CheckCwAssign(cluster_id, size, p[k]))
				{
					// p[k] is the selected path for req., and add p[k] into Path_set.
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[k]));
					return cluster_id;
				}
			}
		}
	}
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
	}
	return ci;
}

int cbrsa::FindCwToAssign_KSP_MF_1(int req_id, int size, int Kpath, PhysicalPath p[])
{
	if(!mf_Cw_seq[size].empty())
	{
		int select_cid = -1;
		int select_k = -1;
		int min_cost = INT_MAX;
		int select_size = INT_MAX;
		int select_hop = INT_MAX;
		for(int k=0; k<Kpath; k++)
		{
			for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
			{
				int cluster_id = it->second;
				if(CheckCwAssign(cluster_id, size, p[k]) && select_size>Cw[cluster_id].size)
				{
					int cost = (Cw[cluster_id].size-size) * (p[k].node.size()-1);
					if(cost < min_cost)
					{
						min_cost = cost;
						select_cid = cluster_id;
						select_k = k;
						select_size = Cw[cluster_id].size;
						select_hop = p[k].node.size()-1;
					}
				}
			}
			if(select_cid != -1 && k<Kpath-1 && select_hop < p[k+1].node.size()-1)
				break;
		}
		if(select_cid != -1)
		{
			if(size*(p[0].node.size()-1) < min_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[select_k]));
			return select_cid;
		}
	}
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
	}
	return ci;
}

/*
int cbrsa::FindCwToAssign_KSP_MF_2(int req_id, int size, int Kpath, PhysicalPath p[])
{
	if(!mf_Cw_seq[size].empty())
	{
		int select_cid = -1;
		int select_k = -1;
		int min_cost = INT_MAX;
		int select_size = INT_MAX;
		int select_hop = INT_MAX;
		for(int k=0; k<Kpath; k++)
		{
			for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
			{
				int cluster_id = it->second;
				if(CheckCwAssign(cluster_id, size, p[k]) && select_size>Cw[cluster_id].size)
				{
					int cost = (Cw[cluster_id].size-size) * (p[k].node.size()-1);
					if(cost < min_cost)
					{
						min_cost = cost;
						select_cid = cluster_id;
						select_k = k;
						select_size = Cw[cluster_id].size;
						select_hop = p[k].node.size()-1;
					}
				}
			}
			if(select_cid != -1 && k<Kpath-1 && select_hop < p[k+1].node.size()-1)
				break;
		}
		if(select_cid != -1)
		{
			if(size*(p[0].node.size()-1) < min_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[select_k]));
			return select_cid;
		}
	}
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
		return ci;
	}
	else  // no assignment in existing clusters and no cluster creation 
	{
		//cout<<"Try to use multiple clusters"<<endl;
		// find multiple clusters to assign
		int min_cw_size = 1;
		int select_cid = -1;
		int min_multi_size = INT_MAX;
		int select_hop = INT_MAX;
		int min_cost = INT_MAX;
		int select_k = -1;
		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[min_cw_size].begin(); it!=mf_Cw_seq[min_cw_size].end(); it++)
		{
			int cluster_id = it->second;
			if(Cw[cluster_id].size >= size)
				break;
			// find countinuous clusters to assign
			int idx = find(sidx_Cw_seq.begin(), sidx_Cw_seq.end(), cluster_id) - sidx_Cw_seq.begin();  // the index of the cluster_id in sidx_Cw_seq[]
			int multi_size = Cw[cluster_id].size;
			int final_idx = idx+1;
			int check_continue = Cw[cluster_id].sidx + Cw[cluster_id].size;
			while(final_idx < sidx_Cw_seq.size())
			{
				if(check_continue != Cw[sidx_Cw_seq[final_idx]].sidx)  // non-continuous
					break;
				multi_size += Cw[sidx_Cw_seq[final_idx]].size;
				if(multi_size >= size)  // find enough length
					break;
				final_idx++;
			}
			if(multi_size < size || min_multi_size < multi_size)
				continue;
			else
			{
				//cout<<"find enough clusters start at Cw["<<sidx_Cw_seq[idx]<<"] with size "<<multi_size<<endl;
				// conbine the weight of Cw in sidx_Cw_seq[idx,m]
				vector<int> combine_weight(et.edge_num, 1);
				for(int i=idx; i<=final_idx; i++)
					for(int j=0; j<et.edge_num; j++)
						combine_weight[j] = (bool)combine_weight[j] & (bool)Cw[sidx_Cw_seq[i]].weight[j];
				// check whether the path can assign into combine_weight
				for(int k=0; k<Kpath; k++)
				{
					if(CheckCwAssign(&combine_weight[0], p[k]))
					{
						int cost = (multi_size-size) * (p[k].node.size()-1);
						if(cost < min_cost)
						{
							min_cost = cost;
							select_cid = cluster_id;
							min_multi_size = multi_size;
							select_k = k;
							select_hop = p[k].node.size()-1;
							select_multi_cw.clear();
							for(int i=idx; i<=final_idx; i++)
								select_multi_cw.push_back(i);
							break;
						}
					}
				}
			}
		}  // end for
		if(select_cid != -1)
		{
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[select_k]));
			// the Cw set is stored in the select_multi_cw[]
			return select_cid;
		}
	}  // end else
	return -1;
}
*/

int cbrsa::FindCwToAssign_AR_SP_1(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	if(!sidx_Cw_seq.empty())
	{
		int select_cid = -1;
		// all pair shortest path
		Graph my_graph("input/topology.txt");
		YenTopKShortestPathsAlg yenAlg(my_graph, my_graph.get_vertex(sour), my_graph.get_vertex(dest));
		int control_flag = true;
		while(yenAlg.has_next())
		{
			BasePath* p = yenAlg.next();
			if(control_flag && p->Weight() > default_p.len)  // try to create a new cluster
			{
				//cout<<p->Weight()<<","<<default_p.len<<"--> try to create a new cluster."<<endl;
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
				control_flag = false;
			}
			PhysicalPath tmp_p;
			tmp_p.len = p->Weight();
			tmp_p.node.reserve(p->length());				
			for(std::vector<BaseVertex*>::const_iterator pos=p->m_vtVertexList.begin(); pos!=p->m_vtVertexList.end();++pos)
				tmp_p.node.push_back((*pos)->getID());
			for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
			{
				int cluster_id = sidx_Cw_seq[i];
				if(Cw[cluster_id].size < size)  // the cluster size < req's size
					continue;
				if(CheckCwAssign(cluster_id, size, tmp_p))
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,tmp_p));
					return cluster_id;
				}
			}
		}
		return -1;
	}

	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
		return ci;
	}
	return -1;
}

int cbrsa::FindCwToAssign_AR_SP(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	if(!sidx_Cw_seq.empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		//select_p.len = pool.sum>=size?default_p.len:INT_MAX;
		select_p.len = INT_MAX;
		select_p.node.clear();
		vector<int> p_node;
		for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
		{
			int cluster_id = sidx_Cw_seq[i];
			if(Cw[cluster_id].size < size)  // the cluster size < req's size
				continue;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1 && p.len<select_p.len)  // existing path and shorter than previous
			{
				// assign cluster_id to select_cid;
				select_cid = cluster_id;
				// assign p to select_p
				select_p.len = p.len;
				p_node.clear();
				p_node.assign(p.node.begin(), p.node.end());			
			}
		}
		if(select_cid != -1)
		{
			select_p.node.assign(p_node.begin(), p_node.end());
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			return select_cid;
		}
	}

	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_SP_2(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// SP
	if(!sidx_Cw_seq.empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		select_p.len = INT_MAX;
		int select_hop = INT_MAX;
		int min_hop = default_p.node.size()-1;
		for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
		{
			int cluster_id = sidx_Cw_seq[i];
			if(Cw[cluster_id].size < size)  // the cluster size < req's size
				continue;
			if(select_cid != -1 && Cw[select_cid].size*select_hop < Cw[cluster_id].size*min_hop)
				continue;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1 && p.node.size()-1<select_hop)  // existing path and shorter than previous
			{
				// assign cluster_id to select_cid;
				select_cid = cluster_id;
				// assign p to select_p
				select_p.len = p.len;
				select_p.node.assign(p.node.begin(), p.node.end());			
				select_hop = p.node.size()-1;
			}
		}
		if(select_cid != -1)
		{
			int create_cost = size*(default_p.node.size()-1);
			int waste_cost = (Cw[select_cid].size-size)*select_hop;
			if(create_cost < waste_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
			}

			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			return select_cid;
		}
	}

	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_SP_4(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// SP + MF w/ Dijkstra
	if(!sidx_Cw_seq.empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		int select_hop = INT_MAX;
		int select_size = INT_MAX;

		for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
		{
			int cluster_id = sidx_Cw_seq[i];
			if(Cw[cluster_id].size < size || Cw[cluster_id].size > select_size)  // the cluster size < req's size
				continue;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1 && (p.node.size()-1<select_hop || (p.node.size()-1==select_hop && Cw[cluster_id].size < select_size) ))
			{
				// assign cluster_id to select_cid;
				select_cid = cluster_id;
				select_size = Cw[cluster_id].size;
				// assign p to select_p
				select_p.len = p.len;
				select_p.node.assign(p.node.begin(), p.node.end());			
				select_hop = p.node.size()-1;
			}
		}

		if(select_cid != -1)
		{
			int create_cost = size*(default_p.node.size()-1);
			int waste_cost = (Cw[select_cid].size-size)*select_hop;
			if(create_cost < waste_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			//Assign_Analysis(size, select_p.node.size()-1);
			return select_cid;
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_SP_6(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// SP + MF w/ Dijkstra with size_blocking solution
	if(!sidx_Cw_seq.empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		int select_hop = INT_MAX;
		int select_size = INT_MAX;

		for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
		{
			int cluster_id = sidx_Cw_seq[i];
			if(Cw[cluster_id].size < size || Cw[cluster_id].size > select_size)  // the cluster size < req's size
				continue;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1 && (p.node.size()-1<select_hop || (p.node.size()-1==select_hop && Cw[cluster_id].size < select_size) ))
			{
				// assign cluster_id to select_cid;
				select_cid = cluster_id;
				select_size = Cw[cluster_id].size;
				// assign p to select_p
				select_p.len = p.len;
				select_p.node.assign(p.node.begin(), p.node.end());			
				select_hop = p.node.size()-1;
			}
		}

		if(select_cid != -1)
		{
			int create_cost = size*(default_p.node.size()-1);
			int waste_cost = (Cw[select_cid].size-size)*select_hop;
			if(create_cost < waste_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					AddIntoWaitingQueue(req_id, size, ci, default_p);
					return ci;
				}		
			}
			// allocate in select_cid
			AddIntoWaitingQueue(req_id, size, select_cid, select_p);
			return select_cid;
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		AddIntoWaitingQueue(req_id, size, ci, default_p);
		return ci;
	}
	else  // no assignment in existing clusters and no cluster creation 
	{
		//cout<<"Try to use multiple clusters"<<endl;

		// find multiple clusters to assign
		int min_cw_size = 1;  // the min size in mf_Cw_seq[] to search
		int select_cid = -1;  // the cluster id of the starting multiple cluster.
		int min_multi_size = INT_MAX;  // the min total size of the super-cluster.
		PhysicalPath select_p;
		int select_hop = INT_MAX;  // the hop of the select_p.
		vector<int> multi_cw;  // store the multiple cw in result.

		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[min_cw_size].begin(); it!=mf_Cw_seq[min_cw_size].end(); it++)
		{
			int cluster_id = it->second;
			if(Cw[cluster_id].size >= size)
				break;
			// find countinuous clusters to assign
			int idx = find(sidx_Cw_seq.begin(), sidx_Cw_seq.end(), cluster_id) - sidx_Cw_seq.begin();  // the index of the cluster_id in sidx_Cw_seq[]
			int multi_size = Cw[cluster_id].size;
			int check_continue = Cw[cluster_id].sidx + Cw[cluster_id].size;
			int final_idx = idx+1;
			while(final_idx < sidx_Cw_seq.size())
			{
				if(check_continue != Cw[sidx_Cw_seq[final_idx]].sidx)  // non-continuous
					break;
				multi_size += Cw[sidx_Cw_seq[final_idx]].size;  // add size
				check_continue = Cw[sidx_Cw_seq[final_idx]].sidx + Cw[sidx_Cw_seq[final_idx]].size;  // update check_continue
				if(multi_size >= size)  // find enough length
					break;
				final_idx++;
			}
			if(multi_size < size || min_multi_size < multi_size)
				continue;
			else
			{
				//cout<<"find enough clusters start at Cw["<<sidx_Cw_seq[idx]<<"] with size "<<multi_size<<endl;

				// conbine the weight of Cw in sidx_Cw_seq[idx,m]
				vector<int> combine_weight(et.edge_num, 1);
				for(int i=idx; i<=final_idx; i++)
				{
					for(int j=0; j<et.edge_num; j++)
					{
						if(combine_weight[j])
							combine_weight[j] = (bool)combine_weight[j] & (bool)Cw[sidx_Cw_seq[i]].weight[j];
					}
				}
				//ShowVector("combine_weight[]", combine_weight);

				// check whether the path can assign into combine_weight
				PhysicalPath p = DijkstraSP(sour, dest, et.node_num, &combine_weight[0]);
				if(p.node.size()>1 && p.node.size()-1<select_hop)
				{
					// assign cluster_id to select_cid;
					select_cid = cluster_id;
					// assign p to select_p
					select_p.len = p.len;
					select_p.node.assign(p.node.begin(), p.node.end());
					select_hop = select_p.node.size()-1;
					// ---
					min_multi_size = multi_size;
					multi_cw.clear();
					for(int i=idx; i<=final_idx; i++)
					{
						multi_cw.push_back(sidx_Cw_seq[i]);  // add cluster_id into select_multi_cw[]
					}
					//ShowVector("multi_cw[]", multi_cw);
				}
			}
		}  // end for

		if(select_cid != -1)
		{
			vector<int> size_set;
			int remain_size = size;
			for(std::vector<int>::size_type i=0; i!=multi_cw.size(); i++)  // set size_set[]
			{
				size_set.push_back(remain_size>Cw[multi_cw[i]].size?Cw[multi_cw[i]].size:remain_size);
				remain_size -= Cw[multi_cw[i]].size;
			}

			AddIntoWaitingQueue(req_id, size, multi_cw, size_set, select_p);
			return select_cid;
		}

	}  // end else

	return -1;
}

int cbrsa::FindCwToAssign_AR_FF(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	if(!sidx_Cw_seq.empty())
	{
		for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size(); i++)
		{
			int cluster_id = sidx_Cw_seq[i];
			if(Cw[cluster_id].size < size)  // the cluster size < req's size
				continue;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1)  // existing path
			{
				wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p));
				return cluster_id;
			}
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_MF(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// MF
	if(!mf_Cw_seq[size].empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
		{
			int cluster_id = it->second;
			if(Cw[cluster_id].size < size)  // the cluster size < req's size
				continue;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1)  // existing path
			{
				select_cid = cluster_id;  // assign cluster_id to select_cid;
				select_p.len = p.len;  // assign p to select_p
				select_p.node.assign(p.node.begin(), p.node.end());	
				//wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p));
				//return cluster_id;
				break;
			}
		}
		if(select_cid != -1)
		{
			int create_cost = size*(default_p.node.size()-1);
			int waste_cost = (Cw[select_cid].size-size)*(select_p.node.size()-1);
			//int assign_cost = Cw[select_cid].size*select_hop;
			//if(create_cost < waste_cost || assign_cost > 2*create_cost)
			if(create_cost < waste_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			//Assign_Analysis(size, select_p.node.size()-1);
			return select_cid;
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_MF_1(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// MF + SP w/ Dijkstra
	if(!mf_Cw_seq[size].empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		int select_hop = INT_MAX;
		int select_size = INT_MAX;

		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
		{
			int cluster_id = it->second;
			if(Cw[cluster_id].size > select_size)
				break;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			if(p.node.size()>1 && p.node.size()-1<select_hop)
			{
				select_cid = cluster_id;  // assign cluster_id to select_cid;
				select_size = Cw[cluster_id].size;
				select_p.len = p.len;  // assign p to select_p
				select_p.node.assign(p.node.begin(), p.node.end());	
				select_hop = p.node.size()-1;
			}
		}

		if(select_cid != -1)
		{
			//int create_cost = size*(default_p.node.size()-1);
			//int waste_cost = (Cw[select_cid].size-size)*select_hop;
			//int assign_cost = Cw[select_cid].size*select_hop;
			//if(create_cost < waste_cost || assign_cost > 2*create_cost)
			if(size*(default_p.node.size()-1) < (Cw[select_cid].size-size)*(select_p.node.size()-1))
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			//Assign_Analysis(size, select_p.node.size()-1);
			return select_cid;
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_MF_2(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// min cost (MC) w/ Dijkstra
	if(!mf_Cw_seq[size].empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		int select_hop = INT_MAX;
		int select_size = INT_MAX;
		int min_cost = INT_MAX;  // cost value
		int min_hop = default_p.node.size()-1;  // the min path hop of the (sour,dest).

		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
		{
			int cluster_id = it->second;
			//if(Cw[cluster_id].size > select_size)
				//break;
			if((Cw[cluster_id].size-size)*min_hop > min_cost)
				break;
			r_count++;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			int cost = (Cw[cluster_id].size-size) * (p.node.size()-1);
			// with the same size of select_size, the shortet path will lower the cost
			if(p.node.size()>1 && cost<min_cost)
			{
				min_cost = cost;
				select_cid = cluster_id;  // assign cluster_id to select_cid;
				select_p.len = p.len;  // assign p to select_p
				select_p.node.assign(p.node.begin(), p.node.end());	
				select_size = Cw[cluster_id].size;
				select_hop = p.node.size()-1;
			}
			if(!cost && select_hop==min_hop)  // best solution
				break;
		}

		if(select_cid != -1)
		{
			int create_cost = size*(default_p.node.size()-1);
			int waste_cost = (Cw[select_cid].size-size)*select_hop;
			int assign_cost = Cw[select_cid].size*select_hop;
			//if(create_cost < waste_cost || assign_cost > 2*create_cost)
			if(size*(default_p.node.size()-1) < min_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			//Assign_Analysis(size, select_p.node.size()-1);
			return select_cid;
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

int cbrsa::FindCwToAssign_AR_MF_6(int req_id, int size, int sour, int dest, PhysicalPath default_p)
{
	// min cost (MC) + SP w/ Dijkstra
	if(!mf_Cw_seq[size].empty())
	{
		int select_cid = -1;
		PhysicalPath select_p;
		int select_hop = INT_MAX;  // the hop of the path select_p
		int min_cost = INT_MAX;  // cost value
		int min_hop = default_p.node.size()-1;  // the min path hop of the (sour,dest).

		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[size].begin(); it!=mf_Cw_seq[size].end(); it++)
		{
			int cluster_id = it->second;
			if((Cw[cluster_id].size-size)*min_hop > min_cost)
				break;  //continue;
			r_count++;
			PhysicalPath p = DijkstraSP(sour, dest, et.node_num, Cw[cluster_id].weight);
			int cost = (Cw[cluster_id].size-size) * (p.node.size()-1);
			//if(p.node.size()>1 && cost<=min_cost && p.node.size()-1 < select_hop)
			if(p.node.size()>1 && (cost<min_cost || (cost==min_cost && p.node.size()-1 < select_hop)))
			{
				min_cost = cost;
				select_cid = cluster_id;  // assign cluster_id to select_cid;
				select_p.len = p.len;  // assign p to select_p
				select_p.node.assign(p.node.begin(), p.node.end());	
				select_hop = p.node.size()-1;
			}
			if(!cost && select_hop==min_hop)  // best solution
				break;
		}

		if(select_cid != -1)
		{
			if(size*(default_p.node.size()-1) < min_cost)
			{
				// try to create a new cluster
				int ci = CreateWCluster_FF(size);
				if(ci != -1)
				{
					wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
					return ci;
				}
			}
			// allocate in select_cid
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,select_p));
			return select_cid;
		}
	}
	
	// try to create a new cluster
	int ci = CreateWCluster_FF(size);
	if(ci != -1)
	{
		wPath_set.insert(std::pair<int,PhysicalPath>(req_id,default_p));
	}
	return ci;
}

// ------------------

PhysicalPath cbrsa::DijkstraSP(int s, int d, int N, int weight[])
{
	//cout<<"Run Dijkstra's Alg"<<endl;
	//for(int i=0; i<et.edge_num; i++)
		//cout<<weight[i]<<" ";
	//cout<<endl;
	PhysicalPath p;
	p.len = 0;
	p.node.clear();
	// initialization
	vector<int> dist(N,INT_MAX/N);
	vector<int> prev(N,-1);
	vector<int> Q(N-1);
	vector<int> found(N,0);
	dist[s] = 0;
	found[s] = 1;
	int j = 0;
	for(int i=0; i<N; i++)
		if(i != s)
			Q[j++] = i;
	//ShowVector("Q[]",Q);

	//start
	int u = s;
	do
	{
		int min_dist = INT_MAX;
		int select = -1;
		for(std::vector<int>::size_type i = 0; i != Q.size(); i++)
		{
			int v = Q[i];
			if(et.eid[u][v]!=-1 && weight[et.eid[u][v]])
			{
				int alt = dist[u] + this->G.graph[u][v];
				if(alt < dist[v])
				{
					dist[v] = alt;
					prev[v] = u;
				}
			}
			// choose
			if(min_dist > dist[v])
			{
				select = v;
				min_dist = dist[v];
			}
		}
		u = select;
		//cout<<u<<" ";
		if(dist[u]>INT_MAX/(2*N))  // no need to search
			break;
		if(u == d)  // arrive at dest
		{
			vector<int> S;
			while(prev[u]!=-1)
			{
				S.push_back(u);
				u = prev[u];
			}
			S.push_back(s);
			p.len = dist[d];
			p.node.reserve(S.size());
			//cout<<"The shortest path from "<<s<<" to "<<d<<" is: ";
			while(!S.empty())
			{
				p.node.push_back(S.back());
				//cout<<S.back()<<" ";
				S.pop_back();
			}
			//cout<<endl;
			return p;
		}
		found[u] = 1;
		Q.erase(std::remove(Q.begin(), Q.end(), u), Q.end());
	} while(!Q.empty());

	/*
	for(int i=0; i<N; i++)
	{
		cout<<"dist["<<i<<"] = "<<dist[i]<<endl;
		cout<<"prev["<<i<<"] = "<<prev[i]<<endl;
	}
	//exit(0);  // del
	*/
	return p;
}

bool cbrsa::CheckCwAssign(int cluster_id, int size, PhysicalPath p)
{
	for(int h=1; h<p.node.size(); h++)  // for each hop of path p[k]
	{
		//link: e(p[k].node[h], p[k].node[h+1])
		int eid = et.eid[p.node[h-1]][p.node[h]];
		//if(Cw[i].weight[eid] >= size)
		if(!Cw[cluster_id].weight[eid])
		{
			// the path is unavailable
			return false;
		}
	}
	return true;
}

bool cbrsa::CheckCwAssign(int weight[], PhysicalPath p)
{
	for(int h=1; h<p.node.size(); h++)  // for each hop of path p[k]
	{
		//link: e(p[k].node[h], p[k].node[h+1])
		int eid = et.eid[p.node[h-1]][p.node[h]];
		if(!weight[eid])
			return false;
	}
	return true;
}

int cbrsa::CreateWCluster_FF(int size)
{
	// find the free block in pool
	int bi = -1;
	for(int i=0; i<pool.num; i++)
	{
		if(pool.block[i][1]-pool.block[i][0]+1 >= size)  // check the block size is bigger than the neccessary size
		{
			bi = i;
			break;
		}
	}
	//cout<<"bi="<<bi<<endl;

	// the index of the new cluster
	int new_ci = -1;
	if(bi != -1)  // there are enough space to create a new cluster
	{
		// create a new cluster
		new_ci = CreateCw(pool.block[bi][0], size);
		/*
		// test - show the mf_Cw_seq
		for(int i=1; i<=15; i++)
		{
			for(vector< pair<int,int> >::iterator it=mf_Cw_seq[i].begin(); it!=mf_Cw_seq[i].end(); it++)
				cout<<"("<< it->first <<","<< it->second <<") ";
			cout<<endl;
		}
		*/

		// update the pool
		if(Cw[new_ci].sidx+size-1 == pool.block[bi][1])
		{
			// the available bloack is filled
			for(int i=bi+1; i<pool.num; i++)
			{
				pool.block[i-1][0] = pool.block[i][0];
				pool.block[i-1][1] = pool.block[i][1];
			}
			pool.num -= 1;
		}
		else
		{
			pool.block[bi][0] += size;
		}
		pool.sum -= size;

		//cout<<"Create Cw["<<new_ci<<"] at ["<<Cw[new_ci].sidx<<","<<Cw[new_ci].sidx+Cw[new_ci].size-1<<"]"<<endl;
		//ShowPool();
	}
	return new_ci;
}

int cbrsa::CreateCw(int sidx, int size)
{
	// get a unused cluster id
	int cid = candidate_wcluster_id.at(candidate_wcluster_id.size()-1);
	candidate_wcluster_id.pop_back();
	// assign info. into new cluster
	Cw[cid].OpNum = 0;
	Cw[cid].OpSet.clear();
	Cw[cid].sidx = sidx;
	Cw[cid].size = size;
	for(int j=0; j<et.edge_num; j++)
		Cw[cid].weight[j] = size;
	// add new id into searching sequence sidx_Cw_seq (for FF and SP), and sort by sidx
	sidx_Cw_seq.push_back(cid);
	// sorting sidx_Cw_seq
	Sort_by_sidx();
	// add new id into mf_Cw_seq[] (searching seq. for MF)
	AddCwIntoMfCwSeq(cid, size);
	// add new id into Cw_size[]
	Cw_size[size].push_back(cid);

	// analysis
	total_Cw_num++;
	sum_Cw_size += Cw[cid].size;

	/*
	// show mf_Cw_seq[] and Cw_size[].
	cout<<"mf_Cw_seq:"<<endl;
	for(int i=1; i<=ReqMaxSize; i++)
	{
		cout<<"["<<i<<"]: ";
		for(vector< pair<int,int> >::iterator it=mf_Cw_seq[i].begin(); it!=mf_Cw_seq[i].end(); it++)
			cout<<"("<<it->first<<","<<it->second<<") ";
		cout<<endl;
	}
	cout<<"Cw_size:"<<endl;
	for(int i=1; i<=ReqMaxSize; i++)
	{
		cout<<"["<<i<<"]: ";
		for(vector<int>::iterator it=Cw_size[i].begin(); it!=Cw_size[i].end(); it++)
			cout<<*it<<" ";
		cout<<endl;
	}
	*/
	return cid;
}

void cbrsa::AddintoPathSet(int req_id, PhysicalPath p)
{
	wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p));
}

void cbrsa::AssignRequest(int req_id)
{
	// get request's info. from waiting_queue[].
	std::map<int, RequestStruct>::iterator map_it = waiting_queue.find(req_id);
	PhysicalPath p = map_it->second.wp;

	// assignment
	vector<int>::iterator size_it = map_it->second.size_set.begin();
	for(vector<int>::iterator it=map_it->second.wcid_set.begin(); it!=map_it->second.wcid_set.end(); it++, size_it++)
	{
		AssignIntoCw(req_id, *size_it, *it, p);
	}

	// analysis
	if(map_it->second.wcid_set.size() > 1)  // super-cluster.
	{
		
		sum_ctn_cw += map_it->second.wcid_set.size();
		total_super_num++;
	}

	// update the req_info[] and waiting_queue[]
	req_info.insert(make_pair(req_id, map_it->second));  // add req into waiting_queue[].
	waiting_queue.erase(map_it);  // remove req from waiting_queue[].
}

void cbrsa::AssignIntoCw(int req_id, int size, int cid, PhysicalPath p)
{
	for(int i=1; i<p.node.size(); i++)
	{
		int eid = et.eid[p.node[i-1]][p.node[i]];
		Cw[cid].weight[eid] = 0;
	}
	Cw[cid].OpNum++;
	Cw[cid].OpSet.push_back(req_id);
	Cw[cid].SizeBucket[size]++;
}

void cbrsa::AssignRequestIntoCw(int req_id, int size, int cid)
{
	// ===== allocate req into cluster =====
	
	// get path from wPath_set according to req_id
	std::map<int, PhysicalPath>::iterator map_it = wPath_set.find(req_id);
	PhysicalPath p = map_it->second;
	for(int i=1; i<p.node.size(); i++)
	{
		int eid = et.eid[p.node[i-1]][p.node[i]];
		Cw[cid].weight[eid] = 0;
	}
	Cw[cid].OpNum++;
	Cw[cid].OpSet.push_back(req_id);
	Cw[cid].SizeBucket[size]++;
	
	// update the req_info
	RequestStruct req;
	req.MultiCw = false;
	req.wcid = cid;
	req.size = size;
	//req.wp = p;
	req.wp.len = p.len;
	req.wp.node.assign(p.node.begin(), p.node.end());
	req_info.insert(std::pair<int,RequestStruct>(req_id,req));
	
	// remove the wPath in wPath_set
	wPath_set.erase(map_it);
}

/*
void cbrsa::AssignRequestIntoCw_1(int req_id, int size, int cid)  // allocate req into cluster
{
	// get path from wPath_set according to req_id
	std::map<int, PhysicalPath>::iterator map_it = wPath_set.find(req_id);
	PhysicalPath p = map_it->second;

	if(select_multi_cw.empty())  // single cluster
	{
		for(int i=1; i<p.node.size(); i++)
		{
			int eid = et.eid[p.node[i-1]][p.node[i]];
			Cw[cid].weight[eid] = 0;
		}
		Cw[cid].OpNum++;
		Cw[cid].OpSet.push_back(req_id);
		Cw[cid].SizeBucket[size]++;
	
		// update the req_info
		RequestStruct req;
		req.MultiCw = false;
		req.wcid = cid;
		req.size = size;
		//req.wp = p;
		req.wp.len = p.len;
		req.wp.node.assign(p.node.begin(), p.node.end());
		req_info.insert(std::pair<int,RequestStruct>(req_id,req));
	}
	else  // multiple clusters
	{
		for(int i=1; i<p.node.size(); i++)
		{
			int eid = et.eid[p.node[i-1]][p.node[i]];
			for(vector<int>::size_type j=0; j!=select_multi_cw.size(); j++)
			{
				int cluster_id = select_multi_cw[j];
				Cw[cluster_id].weight[eid] = 0;
				Cw[cluster_id].OpNum++;
				Cw[cluster_id].OpSet.push_back(req_id);
				Cw[cid].SizeBucket[Cw[cluster_id].size]++;  // tmp
			}
		}

		// update the req_info
		RequestStruct req;
		req.MultiCw = true;
		//req.wcid = cid;
		req.size = size;
		req.wcid_set.assign(select_multi_cw.begin(), select_multi_cw.end());
		//req.wp = p;
		req.wp.len = p.len;
		req.wp.node.assign(p.node.begin(), p.node.end());
		req_info.insert(std::pair<int,RequestStruct>(req_id,req));

		//cout<<"assign req into multiple clusters: ";
		//for(vector<int>::size_type j=0; j!=req.wcid_set.size(); j++)
			//cout<<req.wcid_set[j]<<" ";
		//cout<<endl;

		// analysis
		sum_ctn_cw += select_multi_cw.size();
		total_super_num++;

		// clear the select_multi_cw[]
		select_multi_cw.clear();
	}

	// remove the wPath in wPath_set
	wPath_set.erase(map_it);
}
*/

void cbrsa::ReleaseRequest(int req_id)
{
	// get request's info.
	std::map<int, RequestStruct>::iterator map_it = req_info.find(req_id);
	PhysicalPath p = map_it->second.wp;

	// delete sub-demand from subC
	vector<int>::iterator size_it = map_it->second.size_set.begin();
	for(vector<int>::iterator it=map_it->second.wcid_set.begin(); it!=map_it->second.wcid_set.end(); it++, size_it++)
	{
		ReleaseFromCw(req_id, *size_it, *it, p);
	}

	// update
	req_info.erase(map_it);
}

void cbrsa::ReleaseFromCw(int req_id, int size, int cluster_id, PhysicalPath p)
{
	// recover the edge weight
	for(int i=1; i<p.node.size(); i++)
	{
		int eid = et.eid[p.node[i-1]][p.node[i]];
		Cw[cluster_id].weight[eid] = Cw[cluster_id].size;
	}

	// update the variables in cluster
	Cw[cluster_id].OpNum--;
	Cw[cluster_id].OpSet.erase(std::remove(Cw[cluster_id].OpSet.begin(), Cw[cluster_id].OpSet.end(), req_id), Cw[cluster_id].OpSet.end());
	Cw[cluster_id].SizeBucket[size]--;

	// check whether the cluster still had other members. If not, release the cluster and update the other relative date structure
	if(Cw[cluster_id].OpNum == 0)  // delete the cluster Cw[cid] because of no op in it.
	{
		candidate_wcluster_id.push_back(cluster_id);
		
		//remove Cw[cluster_id] from sidx_Cw_seq[], mf_Cw_seq[], Cw_size[]
		RemoveCwFromSidxCwSeq(cluster_id);
		RemoveCwFromMfCwSeq(cluster_id, Cw[cluster_id].size);
		RemoveCwFromCwSize(cluster_id, Cw[cluster_id].size);
		
		//cout<<"Release ["<<Cw[Cid].sidx<<","<<Cw[Cid].sidx+Cw[Cid].size-1<<"]"<<endl;
		AddFreeBlockIntoPool(Cw[cluster_id].sidx, Cw[cluster_id].sidx+Cw[cluster_id].size-1);  // update pool
		//ShowPool();
		
		// update the variable in Cw[cluster_id]
		Cw[cluster_id].SizeBucket.assign(ReqMaxSize+1,0);
	}
	/*
	else if(size==Cw[Cid].size && !Cw[Cid].SizeBucket[size])  // cluster contract (change the size of Cw[Cid])
	{
		// find the max size of op in cluster
		int max_size = 0, size_i = size-1;
		while(size_i>0)
		{
			if(Cw[Cid].SizeBucket[size_i])
			{
				max_size = size_i;
				break;
			}
			else
				size_i--;
		}

		// update mf_Cw_seq[], Cw_size[]
		for(int i=1; i<=Cw[Cid].size; i++)
			mf_Cw_seq[i].erase(std::remove(mf_Cw_seq[i].begin(), mf_Cw_seq[i].end(), pair<int,int>(Cw[Cid].size-i,Cid)), mf_Cw_seq[i].end());
		for(int i=1; i<=max_size; i++)
		{
			mf_Cw_seq[i].push_back(std::pair<int,int>(max_size-i,Cid));
			std::sort(mf_Cw_seq[i].begin(),mf_Cw_seq[i].end());
		}
		Cw_size[Cw[Cid].size].erase(std::remove(Cw_size[Cw[Cid].size].begin(), Cw_size[Cw[Cid].size].end(), Cid), Cw_size[Cw[Cid].size].end());
		Cw_size[max_size].push_back(Cid);

		// Release the contract part
		AddFreeBlockIntoPool(Cw[Cid].sidx+max_size, Cw[Cid].sidx+Cw[Cid].size-1);
		
		// change the cluster size
			//cout<<"Contract Cw["<<Cid<<"]'s size from "<<Cw[Cid].size<<" to "<<max_size<<endl;
		Cw[Cid].size = max_size;
		for(int i=0; i<et.edge_num; i++)  // update Cw[Cid].weight[]
			if(Cw[Cid].weight[i])
				Cw[Cid].weight[i] = max_size;
		//ShowPool();
	}
	*/
}

void cbrsa::Release_Request(int req_id)
{
	// get the req's information
	std::map<int, RequestStruct>::iterator map_it = req_info.find(req_id);
	int Cid = map_it->second.wcid;  // the cluster index of the req_id assignment.
	int size = map_it->second.size;  // the size of the req_id.
	
	// recover the edge weight
	PhysicalPath p = map_it->second.wp;
	for(int i=1; i<p.node.size(); i++)
	{
		int eid = et.eid[p.node[i-1]][p.node[i]];
		Cw[Cid].weight[eid] = Cw[Cid].size;
	}

	// remove the req from req_info.
	req_info.erase(map_it);

	// update the member in cluster
	Cw[Cid].OpNum--;
	Cw[Cid].OpSet.erase(std::remove(Cw[Cid].OpSet.begin(), Cw[Cid].OpSet.end(), req_id), Cw[Cid].OpSet.end());

	// update the size_bucket
	Cw[Cid].SizeBucket[size]--;

	// check whether the cluster still had other members.
	// if not, release the cluster and update the other relative date structure
	if(Cw[Cid].OpNum == 0)  // delete the cluster Cw[cid] because of no op in it.
	{
		candidate_wcluster_id.push_back(Cid);
		
		//remove Cw[Cid] from sidx_Cw_seq[], mf_Cw_seq[], Cw_size[]
		sidx_Cw_seq.erase(std::remove(sidx_Cw_seq.begin(), sidx_Cw_seq.end(), Cid), sidx_Cw_seq.end());
		for(int i=1; i<=Cw[Cid].size; i++)
			mf_Cw_seq[i].erase(std::remove(mf_Cw_seq[i].begin(), mf_Cw_seq[i].end(), pair<int,int>(Cw[Cid].size-i,Cid)), mf_Cw_seq[i].end());
		Cw_size[Cw[Cid].size].erase(std::remove(Cw_size[Cw[Cid].size].begin(), Cw_size[Cw[Cid].size].end(), Cid), Cw_size[Cw[Cid].size].end());
		
		// update pool
		//cout<<"Release ["<<Cw[Cid].sidx<<","<<Cw[Cid].sidx+Cw[Cid].size-1<<"]"<<endl;
		AddFreeBlockIntoPool(Cw[Cid].sidx, Cw[Cid].sidx+Cw[Cid].size-1);
		//ShowPool();
		
		// update the SizeBucket
		Cw[Cid].SizeBucket.assign(ReqMaxSize+1,0);
	}
	/*
	else if(size==Cw[Cid].size && !Cw[Cid].SizeBucket[size])  // cluster contract (change the size of Cw[Cid])
	{
		// find the max size of op in cluster
		int max_size = 0, size_i = size-1;
		while(size_i>0)
		{
			if(Cw[Cid].SizeBucket[size_i])
			{
				max_size = size_i;
				break;
			}
			else
				size_i--;
		}

		// update mf_Cw_seq[], Cw_size[]
		for(int i=1; i<=Cw[Cid].size; i++)
			mf_Cw_seq[i].erase(std::remove(mf_Cw_seq[i].begin(), mf_Cw_seq[i].end(), pair<int,int>(Cw[Cid].size-i,Cid)), mf_Cw_seq[i].end());
		for(int i=1; i<=max_size; i++)
		{
			mf_Cw_seq[i].push_back(std::pair<int,int>(max_size-i,Cid));
			std::sort(mf_Cw_seq[i].begin(),mf_Cw_seq[i].end());
		}
		Cw_size[Cw[Cid].size].erase(std::remove(Cw_size[Cw[Cid].size].begin(), Cw_size[Cw[Cid].size].end(), Cid), Cw_size[Cw[Cid].size].end());
		Cw_size[max_size].push_back(Cid);

		// Release the contract part
		AddFreeBlockIntoPool(Cw[Cid].sidx+max_size, Cw[Cid].sidx+Cw[Cid].size-1);
		
		// change the cluster size
			//cout<<"Contract Cw["<<Cid<<"]'s size from "<<Cw[Cid].size<<" to "<<max_size<<endl;
		Cw[Cid].size = max_size;
		for(int i=0; i<et.edge_num; i++)  // update Cw[Cid].weight[]
			if(Cw[Cid].weight[i])
				Cw[Cid].weight[i] = max_size;
		//ShowPool();
	}
	*/
}

/*
void cbrsa::Release_Request_1(int req_id)
{
	// get the req's information
	std::map<int, RequestStruct>::iterator map_it = req_info.find(req_id);
	if(map_it->second.MultiCw)  // multiple cluster
	{
		for(vector<int>::size_type idx=0; idx!=map_it->second.wcid_set.size(); idx++)
		{
			int Cid = map_it->second.wcid_set[idx];
			int size = Cw[Cid].size;  // tmp
			// recover the edge weight
			PhysicalPath p = map_it->second.wp;
			for(int i=1; i<p.node.size(); i++)
			{
				int eid = et.eid[p.node[i-1]][p.node[i]];
				Cw[Cid].weight[eid] = Cw[Cid].size;
			}

			// update the member in cluster
			Cw[Cid].OpNum--;
			Cw[Cid].OpSet.erase(std::remove(Cw[Cid].OpSet.begin(), Cw[Cid].OpSet.end(), req_id), Cw[Cid].OpSet.end());

			Cw[Cid].SizeBucket[Cw[Cid].size]--;

			// check whether the cluster still had other members. If not, release the cluster and update the other relative date structure
			if(Cw[Cid].OpNum == 0)  // delete the cluster Cw[cid] because of no op in it.
			{
				candidate_wcluster_id.push_back(Cid);
		
				//remove Cw[Cid] from sidx_Cw_seq[], mf_Cw_seq[], Cw_size[]
				sidx_Cw_seq.erase(std::remove(sidx_Cw_seq.begin(), sidx_Cw_seq.end(), Cid), sidx_Cw_seq.end());
				for(int i=1; i<=Cw[Cid].size; i++)
					mf_Cw_seq[i].erase(std::remove(mf_Cw_seq[i].begin(), mf_Cw_seq[i].end(), pair<int,int>(Cw[Cid].size-i,Cid)), mf_Cw_seq[i].end());
				Cw_size[Cw[Cid].size].erase(std::remove(Cw_size[Cw[Cid].size].begin(), Cw_size[Cw[Cid].size].end(), Cid), Cw_size[Cw[Cid].size].end());
		
				// update pool
				//cout<<"Release ["<<Cw[Cid].sidx<<","<<Cw[Cid].sidx+Cw[Cid].size-1<<"]"<<endl;
				AddFreeBlockIntoPool(Cw[Cid].sidx, Cw[Cid].sidx+Cw[Cid].size-1);
				//ShowPool();
			}
		}	
	}
	else  // single cluster assignment
	{
		int Cid = map_it->second.wcid;  // the cluster index of the req_id assignment.
		int size = map_it->second.size;  // the size of the req_id.
	
		// recover the edge weight
		PhysicalPath p = map_it->second.wp;
		for(int i=1; i<p.node.size(); i++)
		{
			int eid = et.eid[p.node[i-1]][p.node[i]];
			Cw[Cid].weight[eid] = Cw[Cid].size;
		}

		// update the member in cluster
		Cw[Cid].OpNum--;
		Cw[Cid].OpSet.erase(std::remove(Cw[Cid].OpSet.begin(), Cw[Cid].OpSet.end(), req_id), Cw[Cid].OpSet.end());

		Cw[Cid].SizeBucket[Cw[Cid].size]--;

		// check whether the cluster still had other members.
		// if not, release the cluster and update the other relative date structure
		if(Cw[Cid].OpNum == 0)  // delete the cluster Cw[cid] because of no op in it.
		{
			candidate_wcluster_id.push_back(Cid);
		
			//remove Cw[Cid] from sidx_Cw_seq[], mf_Cw_seq[], Cw_size[]
			sidx_Cw_seq.erase(std::remove(sidx_Cw_seq.begin(), sidx_Cw_seq.end(), Cid), sidx_Cw_seq.end());
			for(int i=1; i<=Cw[Cid].size; i++)
				mf_Cw_seq[i].erase(std::remove(mf_Cw_seq[i].begin(), mf_Cw_seq[i].end(), pair<int,int>(Cw[Cid].size-i,Cid)), mf_Cw_seq[i].end());
			Cw_size[Cw[Cid].size].erase(std::remove(Cw_size[Cw[Cid].size].begin(), Cw_size[Cw[Cid].size].end(), Cid), Cw_size[Cw[Cid].size].end());
		
			// update pool
			//cout<<"Release ["<<Cw[Cid].sidx<<","<<Cw[Cid].sidx+Cw[Cid].size-1<<"]"<<endl;
			AddFreeBlockIntoPool(Cw[Cid].sidx, Cw[Cid].sidx+Cw[Cid].size-1);
			//ShowPool();
		}
	}

	// remove the req from req_info.
	req_info.erase(map_it);
}
*/

void cbrsa::Sort_by_sidx()
{
	// sorting the elements in sidx_Cw_seq
	for(std::vector<int>::size_type i = 0; i != sidx_Cw_seq.size()-1; i++)
	{
		for(std::vector<int>::size_type j = i; j != sidx_Cw_seq.size(); j++)
		{
			if(Cw[sidx_Cw_seq[i]].sidx > Cw[sidx_Cw_seq[j]].sidx)
				swap(sidx_Cw_seq[i], sidx_Cw_seq[j]);
		}
	}
}

void cbrsa::ShowCw()
{
	// show cluster's condition
	cout<<"The cluster's condition: "<<endl;
	for(int i=0; i<sidx_Cw_seq.size(); i++)
	{
		int cid = sidx_Cw_seq[i];
		cout<<"Cw["<<cid<<"]: ["<<Cw[cid].sidx<<","<<Cw[cid].sidx+Cw[cid].size-1<<"], size="<<Cw[cid].size<<" num="<<Cw[cid].OpNum<<endl;
		// show OpSet[]
	}
}

// ---------------------------------------------

// ----- Pool Function ----------
void cbrsa::AddFreeBlockIntoPool(int s1, int s2)
{
	pool.sum += s2-s1+1;
	if(pool.num == 0)
	{
		pool.block[0][0] = s1;
		pool.block[0][1] = s2;
		pool.num++;
	}
	else if(s2 < pool.block[0][0])
	{
		if(s2+1 == pool.block[0][0])  // combine
		{
			pool.block[0][0] = s1;
		}
		else  // move the existing elements
		{
			for(int i=pool.num; i>0; i--)
			{
				pool.block[i][0] = pool.block[i-1][0];
				pool.block[i][1] = pool.block[i-1][1];
			}
			pool.num++;
			pool.block[0][0] = s1;
			pool.block[0][1] = s2;
		}
	}
	else if(s1 > pool.block[pool.num-1][1])
	{
		if(s1-1 == pool.block[pool.num-1][1])  // combine
		{
			pool.block[pool.num-1][1] = s2;
		}
		else
		{
			pool.block[pool.num][0] = s1;
			pool.block[pool.num][1] = s2;
			pool.num++;
		}
	}
	else
	{
		for(int i=1; i<pool.num; i++)
		{
			if(s1>pool.block[i-1][1] && s2<pool.block[i][0])
			{
				// add [s1,s2] into block[i]
				if(s1-1==pool.block[i-1][1] && s2+1==pool.block[i][0])
				{
					// combine the left and right
					pool.block[i-1][1] = pool.block[i][1];
					for(int j=i+1; j<pool.num; j++)
					{
						pool.block[j-1][0] = pool.block[j][0];
						pool.block[j-1][1] = pool.block[j][1];
					}
					pool.num--;
				}
				else if(s1-1==pool.block[i-1][1])
				{
					// only combine into the left
					pool.block[i-1][1] = s2;
				}
				else if(s2+1==pool.block[i][0])
				{
					// onlu combine into the right
					pool.block[i][0] = s1;
				}
				else
				{
					// add [s1,s2] into block[i], and move the others
					for(int j=pool.num; j>i; j--)
					{
						pool.block[j][0] = pool.block[j-1][0];
						pool.block[j][1] = pool.block[j-1][1];
					}
					pool.num++;
					pool.block[i][0] = s1;
					pool.block[i][1] = s2;
				}
				break;
			}
		}
	}
}

void cbrsa::ShowPool()
{
	cout<<"The pool condition:"<<endl;
	for(int i=0; i<pool.num; i++)
	{
		cout<<"["<<pool.block[i][0]<<","<<pool.block[i][1]<<"] ";
	}
	cout<<"sum="<<pool.sum;
	cout<<endl;
}

// ------------------------------

// ----- Analysis -----------------

void cbrsa::GetAnalysis_1(cb_Analysis *cba)
{
	cba->r_count += this->r_count;
	cba->size_block += (double)size_block;

	cba->avg_cw_size += (double)sum_Cw_size/(double)total_Cw_num;
	if(total_super_num)
		cba->avg_continuous_cw_assign += (double)sum_ctn_cw/(double)total_super_num;
}

// --------------------------------

// ----- Defragmentation ----------
int cbrsa::DefragProcess(int req_id, int size, int flag, int sour, int dest, int K, PhysicalPath p[])
{
	if(pool.sum >= size)
	{
		// try to move clusters, and create a new cluster
		int ci = DF_MoveAndCreate(size);
		if(ci != -1)
		{
			wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
			return ci;
		}
	}
	else if(pool.sum && pool.sum < size)  // 0 < pool.sum < size
	{
		if(!flag)
		{
			int ci = DF_MoveAndExtend(req_id, size, K, p);  // KSP
			if(ci != -1)
				return ci;
		}
		else
		{
			int ci = DF_MoveAndExtend(req_id, size, sour, dest);  // AR
			if(ci != -1)
				return ci;
		}
	}
	/*
	if(!pool.sum)  // pool.sum == 0
		DF_CwCombination();

	int count = 2;
	while(count)
	{
		if(pool.sum >= size)
		{
			// try to move clusters, and create a new cluster
			int ci = DF_MoveAndCreate(size);
			if(ci != -1)
			{
				wPath_set.insert(std::pair<int,PhysicalPath>(req_id,p[0]));
				return ci;
			}
		}
		else if(pool.sum && pool.sum < size)  // 0 < pool.sum < size
		{
			if(!flag)
			{
				int ci = DF_MoveAndExtend(req_id, size, K, p);  // KSP
				if(ci != -1)
					return ci;
			}
			else
			{
				DF_MoveAndExtend(size, sour, dest);  // AR
			}
		}
		count--;
		if(count)
			DF_CwCombination();
	}
	*/
	return -1;
}

int cbrsa::DF_MoveAndCreate(int size)
{
	//cout<<"\nMoveAndCreate"<<endl;
	// sidx_Cw_seq: the Cw had been sorted out in this vector
	// get b_s
	vector<int> b_s;
	b_s.reserve(sidx_Cw_seq.size()+1);
	b_s.push_back(Cw[sidx_Cw_seq.front()].sidx);
	for(int i=1; i<sidx_Cw_seq.size(); i++)
		b_s.push_back(Cw[sidx_Cw_seq[i]].sidx - (Cw[sidx_Cw_seq[i-1]].sidx+Cw[sidx_Cw_seq[i-1]].size));
	b_s.push_back(Fmax - (Cw[sidx_Cw_seq.back()].sidx+Cw[sidx_Cw_seq.back()].size));
	//ShowVector("b_s", b_s);

	/*
	ShowCw();
	ShowPool();
	*/

	// DF
	int min_cost_idx = -1;
	vector<int> TC(b_s.size(), 0);
	vector<int> usage(b_s.size());
	int LM, RM;
	LM = RM = 0;
	for(int i=0; i<b_s.size(); i++)  // for each block in b_s
	{
		int ES = size-b_s[i];
		int lpi = i-1, rpi = i;
		int left_demands, right_demands;
		int left_moves, right_moves;
		int total_cost = 0;
		left_demands = right_demands =0;
		left_moves = right_moves = 0;
		vector<int> u(b_s.size());

		while(lpi>=0 && rpi<sidx_Cw_seq.size() && ES)
		{
			int llm = min(b_s[lpi], ES);  // local left move
			int left_cost = (Cw[sidx_Cw_seq[lpi]].OpNum + left_demands)*llm;
			int lrm = min(b_s[rpi+1], ES);  // local right move
			int right_cost = (Cw[sidx_Cw_seq[rpi]].OpNum + right_demands)*lrm;
			
			if(left_cost<=right_cost)
			{
				left_demands += Cw[sidx_Cw_seq[lpi]].OpNum;
				left_moves += llm;
				total_cost += left_cost;
				ES -= llm;
				u[lpi] = llm;
				lpi--;
			}
			else
			{
				right_demands += Cw[sidx_Cw_seq[rpi]].OpNum;
				right_moves += lrm;
				total_cost += right_cost;
				ES -= lrm;
				u[rpi+1] = lrm;
				rpi++;
			}
		}  // end while
		
		while(ES && lpi>=0)  // left side
		{
			int llm = min(b_s[lpi], ES);  // local left move
			int left_cost = (Cw[sidx_Cw_seq[lpi]].OpNum + left_demands)*llm;
			left_demands += Cw[sidx_Cw_seq[lpi]].OpNum;
			left_moves += llm;
			total_cost += left_cost;
			ES -= llm;
			u[lpi] = llm;
			lpi--;
		}
		
		while(ES && rpi<sidx_Cw_seq.size())  // right side
		{
			int lrm = min(b_s[rpi+1], ES);  // local right move
			int right_cost = (Cw[sidx_Cw_seq[rpi]].OpNum + right_demands)*lrm;
			right_demands += Cw[sidx_Cw_seq[rpi]].OpNum;
			right_moves += lrm;
			total_cost += right_cost;
			ES -= lrm;
			u[rpi+1] = lrm;
			rpi++;
		}

		// check the min cost
		TC[i] = total_cost;
		if(min_cost_idx==-1 || TC[i]<TC[min_cost_idx])
		{
			min_cost_idx = i;
			usage.clear();
			usage.assign(u.begin(), u.end());
			LM = left_moves;
			RM = right_moves;
		}
	}

	//ShowVector("TC[]", TC);
	//cout<<"min_cost_idx = "<<min_cost_idx<<endl;
	//ShowVector("usage[]", usage);
	//cout<<"LM = "<<LM<<", RM = "<<RM<<endl;

	// moving clusters
	int pivot = min_cost_idx;
	int left_tmp = 0, right_tmp = RM;
	int pre_idx = 0;
	pool.num = 0;

	for(int i=0; i<pivot; i++)
	{
		int cid = sidx_Cw_seq[i];
		left_tmp += usage[i];
		Cw[cid].sidx -= left_tmp;
		if(b_s[i]-usage[i])
		{
			// add the free block into pool
			pool.block[pool.num][0] = pre_idx;
			pool.block[pool.num][1] = Cw[cid].sidx-1;
			pool.num++;
		}
		pre_idx = Cw[cid].sidx + Cw[cid].size;
	}
	
	// the free block pivot would be used for create a new cluster
	int new_Cw_sidx = pre_idx;
	//cout<<"The new cluster would in ["<<new_Cw_sidx<<","<<new_Cw_sidx+size-1<<"]."<<endl;
	if(pivot!=usage.size()-1)
	{
		Cw[sidx_Cw_seq[pivot]].sidx += right_tmp;
		right_tmp -= usage[pivot+1];
		//new_Cw_sidx = pre_idx;
		pre_idx = Cw[sidx_Cw_seq[pivot]].sidx + Cw[sidx_Cw_seq[pivot]].size;
	
		for(int i=pivot+1; i<sidx_Cw_seq.size(); i++)
		{
			int cid = sidx_Cw_seq[i];
			Cw[cid].sidx += right_tmp;
			right_tmp -= usage[i+1];
			if(b_s[i]-usage[i])
			{
				// add the free block into pool
				pool.block[pool.num][0] = pre_idx;
				pool.block[pool.num][1] = Cw[cid].sidx-1;
				pool.num++;
			}
			pre_idx = Cw[cid].sidx + Cw[cid].size;
		}
		if(b_s.back()-usage.back())
		{
			pool.block[pool.num][0] = pre_idx;
			pool.block[pool.num][1] = Fmax-1;
			pool.num++;
		}
	}
	
	// create a new cluster
	int new_ci = CreateCw(new_Cw_sidx, size);
	pool.sum -= size;
	//cout<<"Create a new cluster Cw["<<new_ci<<"] in ["<<new_Cw_sidx<<","<<new_Cw_sidx+size-1<<"]."<<endl;

	// show cluster's condition
	//ShowCw();
	//ShowPool();

	return new_ci;  // return the cluster id of the new created cluster.
}

void cbrsa::DF_CwCombination()
{
	//
	cout<<"\nDF_CwCombination"<<endl;
	// show cluster's condition
	ShowCw();
	ShowPool();
	for(int i=ReqMaxSize; i>0; i--)
	{
		cout<<"Cw_size["<<i;
		ShowVector("]",Cw_size[i]);
	}
	//

	int free_FS = 0;
	int combination_ratio = 1;

	for(int ratio=0; ratio<=combination_ratio; ratio++)
	{
		for(int size_i=ReqMaxSize; size_i>ratio; size_i--)
		{
			int size_j = size_i - ratio;
			//cout<<"size_i="<<size_i<<", size_j="<<size_j<<endl;

			if(Cw_size[size_i].empty()||Cw_size[size_j].empty())
				continue;

			for(vector<int>::iterator i=Cw_size[size_i].begin(); i!=Cw_size[size_i].end();)
			{
				vector<int>::iterator j = (size_i==size_j) ? i+1 : Cw_size[size_j].begin();
				for(j; j!=Cw_size[size_j].end();)
				{
					//cout<<"cid_i="<<*i<<", cid_j="<<*j<<endl;
					if(combineCw(*i,*j))
					{
						int cid_i = *i, cid_j = *j;
						//cout<<"Cw["<<cid_i<<"] can combine with Cw["<<cid_j<<"]"<<endl;
						// add Cw[cid_j] into Cw[cid_i]
						Cw[cid_i].OpNum += Cw[cid_j].OpNum;
						Cw[cid_i].OpSet.insert(Cw[cid_i].OpSet.end(), Cw[cid_j].OpSet.begin(), Cw[cid_j].OpSet.end());
						for(int m=0; m<Cw[cid_j].OpSet.size(); m++)  // change the cluster of the op in the Cw[cid_j].OpSet
							req_info.find(Cw[cid_j].OpSet[m])->second.wcid = cid_i;
						//for(int m=1; m<=ReqMaxSize; m++)
							//Cw[cid_i].SizeBucket[m] += Cw[cid_j].SizeBucket[m];
						for(int m=0; m<et.edge_num; m++)
							Cw[cid_i].weight[m] = Cw[cid_i].weight[m]&Cw[cid_j].weight[m]==0?0:Cw[cid_i].size;

						// update the pool
						AddFreeBlockIntoPool(Cw[cid_j].sidx, Cw[cid_j].sidx+Cw[cid_j].size-1);
						//ShowPool();

						// update free_FS
						free_FS += Cw[cid_j].size;

						// remove cid_j from sidx_Cw_seq[], mf_Cw_seq[]
						RemoveCwFromSidxCwSeq(cid_j);
						RemoveCwFromMfCwSeq(cid_j, Cw[cid_j].size);

						// clear Cw[*j]
						Cw[cid_j].OpNum = 0;
						Cw[cid_j].OpSet.clear();
						//Cw[cid_j].SizeBucket.assign(ReqMaxSize+1,0);
						candidate_wcluster_id.push_back(cid_j);

						// remove Cw[*j] in Cw_size[size_j]
						j = Cw_size[size_j].erase(j);
					}
					else
					{
						++j;
					}
				}
				++i;
			}
		}
	}

	//
	if(free_FS)
	{
		ShowCw();
		ShowPool();
	}
	//
}

int cbrsa::DF_MoveAndExtend(int req_id, int size, int K, PhysicalPath p[])
{
	//cout<<"DF_MoveAndExtend_KSP"<<endl;
	//cout<<"Want "<<size<<" FSs"<<endl;
	int lower = size-pool.sum;  // upper = size
	// if no cluster's size is in the range[lower, size-1], the extension failed
	int extend_flag = false;
	for(int i=lower; i<size; i++)
		if(!Cw_size[i].empty())
		{
			extend_flag = true;
			break;
		}
	if(!extend_flag)
	{
		//cout<<"No cluster can be extended."<<endl;
		return -1;
	}

	// sidx_Cw_seq: the Cw had been sorted out in this vector
	// get b_s
	vector<int> b_s;  // the block seq. which block is between the adjancent Cw[]
	b_s.reserve(sidx_Cw_seq.size()+1);
	vector<int> c_l(Fmax,-1);  // the index of the Cw in sidx_Cw_seq
	b_s.push_back(Cw[sidx_Cw_seq.front()].sidx);
	c_l[sidx_Cw_seq.front()] = 0;
	for(int i=1; i<sidx_Cw_seq.size(); i++)
	{
		b_s.push_back(Cw[sidx_Cw_seq[i]].sidx - (Cw[sidx_Cw_seq[i-1]].sidx+Cw[sidx_Cw_seq[i-1]].size));
		c_l[sidx_Cw_seq[i]] = i;
	}
	b_s.push_back(Fmax - (Cw[sidx_Cw_seq.back()].sidx+Cw[sidx_Cw_seq.back()].size));
	//ShowVector("b_s", b_s);
	//ShowVector("c_l", c_l);
	// show cluster's condition
	//ShowCw();
	//ShowPool();

	int pivot = -1;
	int select_k = -1;  // the accept KSP p[select_k]
	vector<int> usage(b_s.size());
	int move_flag = true;
	int LM = 0, RM = 0;
	// find the usage
	for(int size_i=size-1; size_i>=lower; size_i--)
	{
		if(Cw_size[size_i].empty())
			continue;
		//cout<<"size_i = "<<size_i<<endl;
		//ShowVector("Cw_size[size_i]", Cw_size[size_i]);
		for(vector<int>::iterator it=Cw_size[size_i].begin(); it!=Cw_size[size_i].end(); it++)  // for each possible size of cid
		{
			int cid = *it;
			// check whether there exists path in Cw
			for(int k=0; k<K; k++)
			{
				if(CheckCwAssign(cid,size,p[k]))
				{
					select_k = k;
					break;
				}
			}
			if(select_k == -1)
				continue;
			//cout<<"select p["<<select_k<<"] in Cw["<<cid<<"]."<<endl;
			
			pivot = cid;

			// try to extend Cw[cid]
			int idx = c_l[cid];
			int lpi = idx, rpi = idx+1;
			int ES = size-Cw[cid].size-b_s[idx+1];
			int left_demands, right_demands;
			int left_moves, right_moves;
			int total_cost = 0;
			left_demands = right_demands =0;
			left_moves = right_moves = 0;
			vector<int> u(b_s.size());
			//cout<<"ES="<<ES<<endl;

			if(ES<=0)  // ES<=0, the block is enough for Cw to extend itselt.
			{
				move_flag = false;
			}
			else  // ES>0, need to move other clusters
			{
				move_flag = true;
				u[idx+1] = b_s[idx+1];
				while(lpi>=0 && rpi<sidx_Cw_seq.size() && ES)
				{
					int llm = min(b_s[lpi], ES);  // local left move
					int left_cost = (Cw[sidx_Cw_seq[lpi]].OpNum + left_demands)*llm;
					int lrm = min(b_s[rpi+1], ES);  // local right move
					int right_cost = (Cw[sidx_Cw_seq[rpi]].OpNum + right_demands)*lrm;
			
					if(left_cost<=right_cost)
					{
						left_demands += Cw[sidx_Cw_seq[lpi]].OpNum;
						left_moves += llm;
						total_cost += left_cost;
						ES -= llm;
						u[lpi] = llm;
						lpi--;
					}
					else
					{
						right_demands += Cw[sidx_Cw_seq[rpi]].OpNum;
						right_moves += lrm;
						total_cost += right_cost;
						ES -= lrm;
						u[rpi+1] = lrm;
						rpi++;
					}
				}  // end while
		
				while(ES && lpi>=0)  // left side
				{
					int llm = min(b_s[lpi], ES);  // local left move
					int left_cost = (Cw[sidx_Cw_seq[lpi]].OpNum + left_demands)*llm;
					left_demands += Cw[sidx_Cw_seq[lpi]].OpNum;
					left_moves += llm;
					total_cost += left_cost;
					ES -= llm;
					u[lpi] = llm;
					lpi--;
				}
		
				while(ES && rpi<sidx_Cw_seq.size())  // right side
				{
					int lrm = min(b_s[rpi+1], ES);  // local right move
					int right_cost = (Cw[sidx_Cw_seq[rpi]].OpNum + right_demands)*lrm;
					right_demands += Cw[sidx_Cw_seq[rpi]].OpNum;
					right_moves += lrm;
					total_cost += right_cost;
					ES -= lrm;
					u[rpi+1] = lrm;
					rpi++;
				}
				usage.clear();
				usage.assign(u.begin(), u.end());
				LM = left_moves;
				RM = right_moves;
			}

			break;

		}  // end it

		if(pivot != -1)
			break;

	}  // end size_i

	if(pivot == -1)
	{
		//cout<<"No path can found in extended clusters."<<endl;
		return -1;
	}
	//ShowVector("usage[]", usage);
	//cout<<"LM = "<<LM<<", RM = "<<RM<<endl;

	int pivot_idx = c_l[pivot];
	//cout<<"pivot="<<pivot<<", pivot_idx="<<pivot_idx<<endl;
	//cout<<"move_flag="<<move_flag<<endl;

	if(!move_flag)
	{
		// find the block in pool
		int bi = -1;
		for(int i=0; i<pool.num; i++)
			if(pool.block[i][0] == Cw[pivot].sidx+Cw[pivot].size)
				bi = i;
		//cout<<"bi = "<<bi<<endl;

		// update the pool
		if(Cw[pivot].sidx+size-1 == pool.block[bi][1])
		{
			// the available bloack is filled
			for(int i=bi+1; i<pool.num; i++)
			{
				pool.block[i-1][0] = pool.block[i][0];
				pool.block[i-1][1] = pool.block[i][1];
			}
			pool.num -= 1;
		}
		else
		{
			pool.block[bi][0] += (size-Cw[pivot].size);
		}
		pool.sum -= (size-Cw[pivot].size);
		//cout<<"Cw["<<pivot<<"] would extend from ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+Cw[pivot].size-1<<"] to ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+size-1<<"]."<<endl;
	}
	else
	{
		// moving clusters
		int left_tmp = 0, right_tmp = RM;
		int pre_idx = 0;
		pool.num = 0;  // recount pool

		//cout<<"Cw["<<pivot<<"] would extend from ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+Cw[pivot].size-1<<"] ";

		for(int i=0; i<=pivot_idx; i++)  // move the left side
		{
			int cid = sidx_Cw_seq[i];
			left_tmp += usage[i];
			Cw[cid].sidx -= left_tmp;
			if(b_s[i]-usage[i])
			{
				// add the free block into pool
				pool.block[pool.num][0] = pre_idx;
				pool.block[pool.num][1] = Cw[cid].sidx-1;
				pool.num++;
			}
			pre_idx = Cw[cid].sidx + Cw[cid].size;
		}
		// the extand Cw is sidx_Cw_seq[pivot]
		pre_idx = Cw[pivot].sidx + size;

		//cout<<"to ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+size-1<<"]."<<endl;

		for(int i=pivot_idx+1; i<sidx_Cw_seq.size(); i++)  // move the right side
		{
			int cid = sidx_Cw_seq[i];
			Cw[cid].sidx += right_tmp;
			right_tmp -= usage[i+1];
			if(b_s[i]-usage[i])
			{
				// add the free block into pool
				pool.block[pool.num][0] = pre_idx;
				pool.block[pool.num][1] = Cw[cid].sidx-1;
				pool.num++;
			}
			pre_idx = Cw[cid].sidx + Cw[cid].size;
		}
		if(b_s.back()-usage.back())  // the last block
		{
			pool.block[pool.num][0] = pre_idx;
			pool.block[pool.num][1] = Fmax-1;
			pool.num++;
		}
		pool.sum -= (size-Cw[pivot].size);
	}

	// update the info. of the extended Cw[pivot]
	RemoveCwFromMfCwSeq(pivot, Cw[pivot].size);  // update mf_Cw_seq[]
	AddCwIntoMfCwSeq(pivot, size);
	RemoveCwFromCwSize(pivot, Cw[pivot].size);  // update Cw_size[]
	Cw_size[size].push_back(pivot);

	// extend cluster Cw[pivot_cid]
	Cw[pivot].size = size;
	AddintoPathSet(req_id, p[select_k]);

	// show cluster's condition
	//ShowCw();
	//ShowPool();
	return pivot;
}

int cbrsa::DF_MoveAndExtend(int req_id, int size, int sour, int dest)
{
	//cout<<"DF_MoveAndExtend_SP"<<endl;
	//cout<<"Want "<<size<<" FSs"<<endl;
	int lower = size-pool.sum;  // upper = size
	// if no cluster's size is in the range[lower, size-1], the extension failed
	int extend_flag = false;
	for(int i=lower; i<size; i++)
		if(!Cw_size[i].empty())
		{
			extend_flag = true;
			break;
		}
	if(!extend_flag)
	{
		//cout<<"No cluster can be extended."<<endl;
		return -1;
	}

	// sidx_Cw_seq: the Cw had been sorted out in this vector
	// get b_s
	vector<int> b_s;  // the block seq. which block is between the adjancent Cw[]
	b_s.reserve(sidx_Cw_seq.size()+1);
	vector<int> c_l(Fmax,-1);  // the index of the Cw in sidx_Cw_seq
	b_s.push_back(Cw[sidx_Cw_seq.front()].sidx);
	c_l[sidx_Cw_seq.front()] = 0;
	for(int i=1; i<sidx_Cw_seq.size(); i++)
	{
		b_s.push_back(Cw[sidx_Cw_seq[i]].sidx - (Cw[sidx_Cw_seq[i-1]].sidx+Cw[sidx_Cw_seq[i-1]].size));
		c_l[sidx_Cw_seq[i]] = i;
	}
	b_s.push_back(Fmax - (Cw[sidx_Cw_seq.back()].sidx+Cw[sidx_Cw_seq.back()].size));
	//ShowVector("b_s", b_s);
	//ShowVector("c_l", c_l);
	// show cluster's condition
	//ShowCw();
	//ShowPool();

	int pivot = -1;
	PhysicalPath select_p;
	vector<int> usage(b_s.size());
	int move_flag = true;
	int LM = 0, RM = 0;
	// find the usage
	for(int size_i=size-1; size_i>=lower; size_i--)
	{
		if(Cw_size[size_i].empty())
			continue;
		//cout<<"size_i = "<<size_i<<endl;
		//ShowVector("Cw_size[size_i]", Cw_size[size_i]);
		
		int select_cid = -1;
		select_p.len = INT_MAX;
		// check whether there exists path in Cw, and choose the SP
		for(vector<int>::iterator it=Cw_size[size_i].begin(); it!=Cw_size[size_i].end(); it++)  // for each possible size of cid
		{
			int cid = *it;
			PhysicalPath p = DijkstraSP(sour,dest,et.node_num,Cw[cid].weight);
			if(p.node.size()>1 && p.len<select_p.len)  // existing path and shorter than previous
			{
				// assign cluster_id to select_cid;
				select_cid = cid;
				// assign p to select_p
				select_p.len = p.len;
				select_p.node.clear();
				select_p.node.assign(p.node.begin(), p.node.end());			
			}
			//cout<<"select Cw["<<select_cid<<"]."<<endl;
		}
		if(select_cid == -1)
			continue;

		int cid = select_cid;
		pivot = cid;

		// try to extend Cw[cid]
		int idx = c_l[cid];
		int lpi = idx, rpi = idx+1;
		int ES = size-Cw[cid].size-b_s[idx+1];
		int left_demands, right_demands;
		int left_moves, right_moves;
		int total_cost = 0;
		left_demands = right_demands =0;
		left_moves = right_moves = 0;
		vector<int> u(b_s.size());
		//cout<<"ES="<<ES<<endl;

		if(ES<=0)  // ES<=0, the block is enough for Cw to extend itselt.
		{
			move_flag = false;
		}
		else  // ES>0, need to move other clusters
		{
			move_flag = true;
			u[idx+1] = b_s[idx+1];
			while(lpi>=0 && rpi<sidx_Cw_seq.size() && ES)
			{
				int llm = min(b_s[lpi], ES);  // local left move
				int left_cost = (Cw[sidx_Cw_seq[lpi]].OpNum + left_demands)*llm;
				int lrm = min(b_s[rpi+1], ES);  // local right move
				int right_cost = (Cw[sidx_Cw_seq[rpi]].OpNum + right_demands)*lrm;
			
				if(left_cost<=right_cost)
				{
					left_demands += Cw[sidx_Cw_seq[lpi]].OpNum;
					left_moves += llm;
					total_cost += left_cost;
					ES -= llm;
					u[lpi] = llm;
					lpi--;
				}
				else
				{
					right_demands += Cw[sidx_Cw_seq[rpi]].OpNum;
					right_moves += lrm;
					total_cost += right_cost;
					ES -= lrm;
					u[rpi+1] = lrm;
					rpi++;
				}
			}  // end while
		
			while(ES && lpi>=0)  // left side
			{
				int llm = min(b_s[lpi], ES);  // local left move
				int left_cost = (Cw[sidx_Cw_seq[lpi]].OpNum + left_demands)*llm;
				left_demands += Cw[sidx_Cw_seq[lpi]].OpNum;
				left_moves += llm;
				total_cost += left_cost;
				ES -= llm;
				u[lpi] = llm;
				lpi--;
			}
		
			while(ES && rpi<sidx_Cw_seq.size())  // right side
			{
				int lrm = min(b_s[rpi+1], ES);  // local right move
				int right_cost = (Cw[sidx_Cw_seq[rpi]].OpNum + right_demands)*lrm;
				right_demands += Cw[sidx_Cw_seq[rpi]].OpNum;
				right_moves += lrm;
				total_cost += right_cost;
				ES -= lrm;
				u[rpi+1] = lrm;
				rpi++;
			}
			usage.clear();
			usage.assign(u.begin(), u.end());
			LM = left_moves;
			RM = right_moves;
		}

		if(pivot != -1)
			break;

	}  // end size_i

	if(pivot == -1)
	{
		//cout<<"No path can found in extended clusters."<<endl;
		return -1;
	}
	//ShowVector("usage[]", usage);
	//cout<<"LM = "<<LM<<", RM = "<<RM<<endl;

	int pivot_idx = c_l[pivot];
	//cout<<"pivot="<<pivot<<", pivot_idx="<<pivot_idx<<endl;
	//cout<<"move_flag="<<move_flag<<endl;

	if(!move_flag)
	{
		// find the block in pool
		int bi = -1;
		for(int i=0; i<pool.num; i++)
			if(pool.block[i][0] == Cw[pivot].sidx+Cw[pivot].size)
				bi = i;
		//cout<<"bi = "<<bi<<endl;

		// update the pool
		if(Cw[pivot].sidx+size-1 == pool.block[bi][1])
		{
			// the available bloack is filled
			for(int i=bi+1; i<pool.num; i++)
			{
				pool.block[i-1][0] = pool.block[i][0];
				pool.block[i-1][1] = pool.block[i][1];
			}
			pool.num -= 1;
		}
		else
		{
			pool.block[bi][0] += (size-Cw[pivot].size);
		}
		pool.sum -= (size-Cw[pivot].size);
		//cout<<"Cw["<<pivot<<"] would extend from ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+Cw[pivot].size-1<<"] to ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+size-1<<"]."<<endl;
	}
	else
	{
		// moving clusters
		int left_tmp = 0, right_tmp = RM;
		int pre_idx = 0;
		pool.num = 0;  // recount pool

		//cout<<"Cw["<<pivot<<"] would extend from ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+Cw[pivot].size-1<<"] ";

		for(int i=0; i<=pivot_idx; i++)  // move the left side
		{
			int cid = sidx_Cw_seq[i];
			left_tmp += usage[i];
			Cw[cid].sidx -= left_tmp;
			if(b_s[i]-usage[i])
			{
				// add the free block into pool
				pool.block[pool.num][0] = pre_idx;
				pool.block[pool.num][1] = Cw[cid].sidx-1;
				pool.num++;
			}
			pre_idx = Cw[cid].sidx + Cw[cid].size;
		}
		// the extand Cw is sidx_Cw_seq[pivot]
		pre_idx = Cw[pivot].sidx + size;

		//cout<<"to ["<<Cw[pivot].sidx<<","<<Cw[pivot].sidx+size-1<<"]."<<endl;

		for(int i=pivot_idx+1; i<sidx_Cw_seq.size(); i++)  // move the right side
		{
			int cid = sidx_Cw_seq[i];
			Cw[cid].sidx += right_tmp;
			right_tmp -= usage[i+1];
			if(b_s[i]-usage[i])
			{
				// add the free block into pool
				pool.block[pool.num][0] = pre_idx;
				pool.block[pool.num][1] = Cw[cid].sidx-1;
				pool.num++;
			}
			pre_idx = Cw[cid].sidx + Cw[cid].size;
		}
		if(b_s.back()-usage.back())  // the last block
		{
			pool.block[pool.num][0] = pre_idx;
			pool.block[pool.num][1] = Fmax-1;
			pool.num++;
		}
		pool.sum -= (size-Cw[pivot].size);
	}

	// update the info. of the extended Cw[pivot]
	RemoveCwFromMfCwSeq(pivot, Cw[pivot].size);  // update mf_Cw_seq[]
	AddCwIntoMfCwSeq(pivot, size);
	RemoveCwFromCwSize(pivot, Cw[pivot].size);  // update Cw_size[]
	Cw_size[size].push_back(pivot);

	// extend cluster Cw[pivot_cid]
	Cw[pivot].size = size;
	AddintoPathSet(req_id, select_p);

	// show cluster's condition
	//ShowCw();
	//ShowPool();

	return pivot;
}

bool cbrsa::combineCw(int m, int n)
{
	// check whether Cw[m] can be add into Cw[n]
	for(int i=0; i<et.edge_num; i++)
		if(!Cw[m].weight[i] && !Cw[n].weight[i])
			return false;
	return true;
}

// --------------------------------

// ----- other function ------------------------
void ShowTwoDimArray(int row, int col, int** arr)
{
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
			cout<<arr[i][j]<<" ";
		cout<<endl;
	}
}

void OutputLinkTopoFile(Topology G, string filename)
{
	string filepath = ".\\input\\" + filename;
	ofstream ofs;
	ofs.open(filepath.c_str(), std::ios_base::trunc);
	ofs<<G.node_num<<endl;
	for(int i=0; i<G.node_num; i++)
	{
		for(int j=0; j<G.node_num; j++)
		{
			if(G.graph[i][j])
			{
				ofs<<i<<'\t'<<j<<'\t'<<G.graph[i][j]<<endl;
			}
		}
	}
	ofs.close();
}

int** CreateTwoDimArray(int row, int col, int value)
{
	int** arr;
	arr = new int*[row];
	for(int i=0; i<row; i++)
	{
		arr[i] = new int[col];
		for(int j=0; j<col; j++)
			arr[i][j] = value;
	}
	return arr;
}

EdgeTransformation SetET(Topology G)
{
	EdgeTransformation et;
	et.node_num = G.node_num;
	et.eid = CreateTwoDimArray(et.node_num, et.node_num, -1);
	int N = et.node_num;
	int E = 0;
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			if(G.graph[i][j])
			{
				et.eid[i][j] = E++;
			}
	et.edge_num = E;
	//ShowTwoDimArray(et.node_num,et.node_num,et.eid);
	return et;
}

//template <class T>
void ShowVector(string str, vector<int> vec)
{
	//typename vector<int>::iterator it;
	vector<int>::iterator it;
	cout<<str<<": ";
	for(it=vec.begin(); it<vec.end(); it++)
		cout<<*it<<" ";
	cout<<endl;
}

// ---------------------------------------------
