#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

struct EvalValue
{
	double BBP;
	double BP;
	int MS;  // maximum index
	double usedFS;  // the total used FS in spectrum
	double duration;  // the time durration of the method
	double path_hop;  // the average path hop for each eop
	// -------------
	double r_time;
	double sa_time;
	vector<double> bp_size;
};

class evaluation
{
private:
	int ReqMaxSize;
	int blockOp;  // the amount of blocking op
	int blockBandwidth;  // the amount of blocking bandwidth
	int MS;  // maximum index
	int totalOp;  // the amount of op for evaluation
	int totalBw;  // the amount of Bw
	int usedFS;  // the total used FS in spectrum
	int total_path_hop;
	// ----------
	vector<int> blockEachSize;
	vector<int> totalEachSize;

public:
	evaluation(int);
	void AcceptEvent(int, int, int);
	void BlockEvent(int);
	void GetAnalsis(EvalValue*);
};