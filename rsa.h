#include <iostream>
#include <vector>
#include <map>

using namespace std;

struct Topology
{
	int node_num;
	vector< vector<int> > graph;
};

struct PhysicalPath
{
    int len;
	vector<int> node;
};

struct RequestStruct
{
	int size;
	bool MultiCw;
	int wcid;  // single cluster assignment
	vector<int> wcid_set;  // multiple cluster assignment: store the multi clusters which req assigned in.
	vector<int> size_set;  // store the sub-size of the req.
	PhysicalPath wp;
};

struct Cluster
{
	int sidx;  // the starting index of the cluster
	int size;  // the size of the cluster
	int OpNum;  // the # of Op in cluster
	vector<int> OpSet;  // the set of Op
	int *weight;  // topology in cluster
	vector<int> SizeBucket;  // record the size number of all the op in this cluster
};

struct Pool
{
	int sum;
	int num;
	int **block;  // block[Fmax][2]
	// block[i] = {m,n}: the i-th free block is in spectrum index [m,n].
};

struct EdgeTransformation
{
	// a struct to trandform edge(N1, N2) into edge id eid[N1][N2].
	int node_num;
	int edge_num;
	int** eid;  // the id of the edge e(n1,n2) is eid[n1][n2];
};

struct cb_Analysis
{
	double r_count;
	double size_block;
	// --------
	double avg_cw_size;
	double avg_continuous_cw_assign;
};

class Routing
{
private:
	int Kpath;
	int Node_num;
    //PhysicalPath*** path;  // path[N][N][K]: record each K-paths of each (s,d) pair in path[][][].
	vector<PhysicalPath>** path;
    int **pairK;  // pairK[N][N]: The (s,d) pair have pairK[s][d] paths, pairK[s][d]<=K.
	// ----------
	//vector<PhysicalPath> **path_set;
	//double **path_time;  // record the counting time of the all path in (s,d) pair. (unit:ms)
	//vector<double> **tmp;

public:
    Routing();
    ~Routing();
	int Get_K_value();
	int Get_Path_Hop(PhysicalPath);
	//void AllShortestPath(int,Topology);  // to calcalate all the path and store in pathset.
	void KShortestPath(int, Topology);  // to calculate and store K-SP.
	PhysicalPath* ReturnKPath(int, int);  // return all the k-path of (s,d) pair.
	PhysicalPath ReturnSinglePath(int, int, int);  // return the k-th path of (s,d) pair.
	int GetPathNum(int, int);
	void ShowPath(PhysicalPath);  // show the input path
};

class cbrsa
{
private:
	int Fmax;
	int ReqMaxSize;
	Topology G;  // used for Dijkstra' SP
	Pool pool;  // record the free blocks
	Cluster* Cw;  // the set of working clusters
	vector<int> sidx_Cw_seq;  // teh searching sequence for working cluster.
	vector<int> candidate_wcluster_id;  // the set of working cluster's id.
	vector< std::pair<int,int> >* mf_Cw_seq;  // store the (mf_value,cid) for searching purpose.
	vector<int>* Cw_size;  // store the size of Cw[i].
	map<int, RequestStruct> waiting_queue;  // store the req waiting for assignment.
	map<int, RequestStruct> req_info;  // record the information of the req.
	map<int, PhysicalPath> wPath_set;  // to store the working path of req_id
	EdgeTransformation et;
	// ----
	//vector<int> select_multi_cw;  // used for multiple clusters, stored the cluster_id of multiple cluters

	// ----- analysis ----------
	int r_count;  // the # of times for searching path.
	int size_block;  // the # of blocking because of the size.
	int total_Cw_num;  // the # of the creatd cluster in the whole instance
	int sum_Cw_size;  // the sum of the cluster's size of the created cluster
	int sum_ctn_cw;  // the sum of the # of the continuous clusters used for one req.
	int total_super_num;  // the # of the req. using super_cluster.

	// -------------------------
	
public:
	cbrsa(Topology, int, int);  // constructor
	~cbrsa();  // destructor

	// ----- get/set value -----
	int found_times;
	int Get_Cw_sidx(int);
	PhysicalPath GetwPath(int);
	void AddIntoWaitingQueue(int, int, int, PhysicalPath);
	void AddIntoWaitingQueue(int, int, vector<int>, vector<int>, PhysicalPath);
	// ---------------------

	// ----- other function -----
	void Sort_by_sidx();
	PhysicalPath DijkstraSP(int, int, int, int[]);
	void ShowCw();
	void RemoveCwFromSidxCwSeq(int);
	void AddCwIntoMfCwSeq(int, int);
	void RemoveCwFromMfCwSeq(int, int);
	void RemoveCwFromCwSize(int, int);

	// --------------------------

	// ----- for request -----
	void AssignRequest(int);  // the interface for req in main.cpp to assign into cluster(s).
	void AssignIntoCw(int, int, int, PhysicalPath);  // assign (sub-)demand into single cluster.
	void ReleaseRequest(int);  // the interface for req in main.cpp to be dropped from cluster(s).
	void ReleaseFromCw(int, int, int, PhysicalPath);  // drop (sub-)demand from single cluster.
	// -----------------------

	// ----- SA method -----
	int FindCwToAssign_KSP_SP(int, int, int, PhysicalPath[]);
	int FindCwToAssign_KSP_FF(int, int, int, PhysicalPath[]);
	int FindCwToAssign_KSP_MF(int, int, int, PhysicalPath[]);
	int FindCwToAssign_KSP_MF_1(int, int, int, PhysicalPath[]);
	//int FindCwToAssign_KSP_MF_2(int, int, int, PhysicalPath[]);
	int FindCwToAssign_AR_SP(int, int, int, int, PhysicalPath);
	int FindCwToAssign_AR_SP_1(int, int, int, int, PhysicalPath);
	int FindCwToAssign_AR_SP_2(int, int, int, int, PhysicalPath);
	int FindCwToAssign_AR_SP_4(int, int, int, int, PhysicalPath);  // SP_MC
	int FindCwToAssign_AR_SP_6(int, int, int, int, PhysicalPath);  // SP_MF_ss
	int FindCwToAssign_AR_FF(int, int, int, int, PhysicalPath);
	int FindCwToAssign_AR_MF(int, int, int, int, PhysicalPath);  // MF
	int FindCwToAssign_AR_MF_1(int, int, int, int, PhysicalPath);
	int FindCwToAssign_AR_MF_2(int, int, int, int, PhysicalPath);
	int FindCwToAssign_AR_MF_6(int, int, int, int, PhysicalPath);  // MC_SP
	bool CheckCwAssign(int, int, PhysicalPath);
	bool CheckCwAssign(int[], PhysicalPath);
	int CreateWCluster_FF(int);
	int CreateCw(int, int);
	void AddintoPathSet(int, PhysicalPath);
	void AssignRequestIntoCw(int, int, int);
	//void AssignRequestIntoCw_1(int, int, int);
	void Release_Request(int);
	//void Release_Request_1(int);

	// ---------------------
	
	// ----- Pool -----
	void AddFreeBlockIntoPool(int, int);
	void ShowPool();
	// ----------------

	// ----- Analysis -----
	void GetAnalysis_1(cb_Analysis*);
	// -------------------

	// ----- DF -----
	int DefragProcess(int, int, int, int, int, int, PhysicalPath[]);
	int DF_MoveAndCreate(int);
	int DF_MoveAndExtend(int,int,int,PhysicalPath[]);  // KSP
	int DF_MoveAndExtend(int,int,int,int);  // AR
	void DF_CwCombination();
	bool combineCw(int,int);
	// void DF_CwRearrange();
	// --------------
};