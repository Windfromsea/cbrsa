#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "evaluation.h"

using namespace std;

evaluation::evaluation(int ReqMaxSize)
{
	this->ReqMaxSize = ReqMaxSize;
	this->blockBandwidth = 0;
	this->blockOp = 0;
	this->MS = 0;
	this->totalBw = 0;
	this->totalOp = 0;
	this->usedFS = 0;
	this->total_path_hop = 0;
	blockEachSize.assign(ReqMaxSize+1, 0);
	totalEachSize.assign(ReqMaxSize+1, 0);
}

void evaluation::AcceptEvent(int size, int sidx, int path_hop)
{
	totalBw += size;
	totalOp++;
	MS = max(MS,sidx+size-1);
	usedFS += size*path_hop;
	total_path_hop += path_hop;
	totalEachSize[size]++;
}

void evaluation::BlockEvent(int size)
{
	blockBandwidth += size;
	blockOp++;
	blockEachSize[size]++;
	totalBw += size;
	totalOp++;
	totalEachSize[size]++;
}

void evaluation::GetAnalsis(EvalValue* e)
{
	e->BBP += (double)blockBandwidth/totalBw;
	e->BP += (double)blockOp/totalOp;
	e->MS += MS;
	e->usedFS += (double)usedFS/(totalOp-blockOp);
	e->path_hop += (double)total_path_hop/(totalOp-blockOp);  // sun of paht_hop / the # of successful op.
	for(int i=1; i<ReqMaxSize+1; i++)
	{
		if(totalEachSize[i])
			e->bp_size[i] += (double)blockEachSize[i]/(double)totalEachSize[i];
	}
	//cout<<blockBandwidth<<" "<<totalBw<<" "<<blockOp<<" "<<totalOp<<" "<<MS<<" "<<usedFS<<" "<<total_path_hop<<endl;
}


