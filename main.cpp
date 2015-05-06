#include "queuing.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <time.h>
// ----- QueryPerformanceCounter -----
#include <windows.h>
#include <Mmsystem.h>
#pragma comment(lib, "Winmm.lib")
// -----------------------------------
//#ifndef RSA_H_
//#define RSA_H_
	#include "rsa.h"
//#endif
#include "evaluation.h"
// -----------------------------------

using namespace std;

// ----- gloabl variable ----------
Routing route;  // used for finding out KSP
double rsa_time;  // to tmp. calculate the duration of rsa
double r_time;
double sa_time;
// --------------------------------
// ----- self definition ----------
enum Action{GENERATE, PROCESS, FINISH, BLOCK};  // req's status when queuing
enum PreRouting{KSP, AR, LastRS};  // the routing strategy: K-SP or Adaptive Routing
enum Cluster_Assignment{FF, SP, MF, Uvalue, LastCAS};  // cluster assignment strategy
//enum Defragment_Strategy{};
PreRouting RS;  // routing strategy
Cluster_Assignment CAS;  // cluster assignment strategy
int DF_flag = 1;
// --------------------------------
// ----- function ----------
void CBRSA(Topology, int, int, int, int);
void Set_parameter(string*, int*, int*, int*, int*, int*, int*, int*, int*);
//void Set_file();
void Generate_request(int, int, int, int);
bool EOP_Setup(cbrsa*, TOKEN*);
void EOP_Release(cbrsa*, TOKEN);
string GetMethodOuptputFile();
void initEvalvalue(EvalValue*, int);
void initCbAalysis(cb_Analysis*, int, int);
Topology Get_topology(string);
void Set_method(int, char*[]);
// -------------------------
// ----- Method ----------
EvalValue ev;  // to store the simulation result of each erlang
cb_Analysis cba;
// -----------------------

int main(int argc, char* argv[])
{
	// ----- variable ----------
	string TopologyFile = "";
	int K = 0;
    int Fmax = 0;
    int TotalReqNum = 0;
    int ReqMaxSize = 0;
    int IniErlang = 0;
    int IntervalErlang = 0;
	int ErlangTimes = 0;
    int EvalTimes = 0;
	// -------------------------

	// set parameter
    Set_parameter(&TopologyFile, &K, &Fmax, &TotalReqNum, &ReqMaxSize, &IniErlang, &IntervalErlang, &ErlangTimes, &EvalTimes);
	// set method
	Set_method(argc, argv);
	// get topology
	Topology mainG = Get_topology(TopologyFile);
	// preditermined routing
	route.KShortestPath(K, mainG);
	/*
	//test for showing the topology.
	for(int i=0; i<mainG.node_num; i++)
		for(int j=0; j<mainG.node_num; j++)
			for(int k=0; k<K; k++)
				route.ShowSingleKPath(i,j,k);
	*/

	// time seed
	srand(time(0));  // srand(time(0)+getpid());
	// simpack initalization
	init_simpack(LINKED);

	string outputfilepath = ".\\analysis\\" + GetMethodOuptputFile();  // output file name
	remove((outputfilepath+".txt").c_str());

	for(int et=0; et<ErlangTimes; et++)
	{
		int current_erlang = IniErlang+IntervalErlang*et;
		initEvalvalue(&ev, ReqMaxSize);  // initial simulation record
		initCbAalysis(&cba, ReqMaxSize, mainG.node_num);
		cout<<"----- Erlang = "<<current_erlang<<" -----"<<endl;
		
		for(int instance=0; instance<EvalTimes; instance++)
		{
			cout<<"Intance "<<instance<<endl;
			rsa_time = r_time = sa_time = 0;

			CBRSA(mainG, Fmax, TotalReqNum, ReqMaxSize, current_erlang);
			
			ev.duration += rsa_time/TotalReqNum;
			ev.r_time += r_time/TotalReqNum;
			ev.sa_time += sa_time/TotalReqNum;
		}

		// show result in the window
		cout<<"Load: "<<current_erlang<<"----------"<<endl;
		cout<<"BBP = "<<(double)ev.BBP/EvalTimes<<endl;
		cout<<"BP = "<<(double)ev.BP/EvalTimes<<endl;
		cout<<"MS = "<<ev.MS/EvalTimes<<endl;
		cout<<"avg_path_hop = "<<(double)ev.path_hop/EvalTimes<<endl;
		cout<<"avg_FS = "<<(double)ev.usedFS/EvalTimes<<endl;
		printf("avg_rsa_time = %llf ms\n", ev.duration/EvalTimes);
		printf("avg_r_time = %llf ms\n", ev.r_time/EvalTimes);
		printf("avg_sa_time = %llf ms\n", ev.sa_time/EvalTimes);
		cout<<"r_count = "<<cba.r_count/EvalTimes<<endl;

		// Output Analsis Result
		ofstream ofs;
		ofs.open((outputfilepath+".txt").c_str(), std::ios_base::app);
		// output all the measure in the output file
		ofs<<current_erlang<<'\t';
		ofs<<(double)ev.BBP/EvalTimes<<'\t';
		ofs<<(double)ev.BP/EvalTimes<<'\t';
		ofs<<ev.MS/EvalTimes<<'\t';
		ofs<<(double)ev.usedFS/EvalTimes<<'\t';
		ofs<<(double)ev.path_hop/EvalTimes<<'\t';
		ofs<<(double)ev.duration/EvalTimes<<'\t';
		ofs<<(double)ev.r_time/EvalTimes<<'\t';
		ofs<<(double)ev.sa_time/EvalTimes<<endl;
		ofs.close();

		/*
		ofs.open((outputfilepath+"_size.txt").c_str(), std::ios_base::app);
		ofs<<current_erlang<<'\t';
		for(int i=1; i<ReqMaxSize+1; i++)
			ofs<<ev.bp_size[i]/EvalTimes<<'\t';
		ofs<<endl;
		ofs.close();
		*/

		ofs.open((outputfilepath+"_cba.txt").c_str(), std::ios_base::app);
		ofs<<current_erlang<<'\t';
		ofs<<cba.avg_continuous_cw_assign/EvalTimes<<'\t';
		ofs<<cba.avg_cw_size/EvalTimes<<endl;
		ofs.close();
	}

	cout<<"Hello world!!!"<<endl;
	return 0;
}

void CBRSA(Topology G, int Fmax,int TotalReqNum, int ReqMaxSize, int current_erlang)
{
	TOKEN token;
	int event;
	int current_req_num = 0;
	bool show_detail = false;
	//bool show_detail = true;
	evaluation eval(ReqMaxSize);

	cbrsa cb(G, Fmax, ReqMaxSize);

	// generate request
	Generate_request(G.node_num, ReqMaxSize, current_req_num, current_erlang);
	current_req_num++;

	// couting time (QueryPerformanceCounter)
	LARGE_INTEGER  large_interger;
	QueryPerformanceFrequency(&large_interger);
	double dff = large_interger.QuadPart;
	__int64 tp1;

	while(current_req_num <= TotalReqNum)
	{
		// get the statu of the event
		int event_statu = next_event(&event,&token);

		if(event_statu != ES_OK)
		{
			break;
		}
		switch(event)
		{
			case PROCESS:
			{
				//double timepoint = clock();  // start to count rsa_time of this instance
				QueryPerformanceCounter(&large_interger);
				__int64 tp1 = large_interger.QuadPart;

				// try to establish a new op, and get the rsa result
				bool accept_eop = EOP_Setup(&cb, &token);

				//rsa_time += (double)(clock()-timepoint)/CLOCKS_PER_SEC;  // end to count rsa_time				
				QueryPerformanceCounter(&large_interger);
				tp1 = large_interger.QuadPart - tp1;
				rsa_time += tp1*1000/dff;

				if(accept_eop)
				{
					// add event into schedule
					schedule(FINISH,expntl(1.0),token);

					if(show_detail)
					{
						// test for output the info. of req. assignment
						cout<<"req "<<token.attr[0]<<" is assign."<<endl;
						cout<<"("<<token.attr[1]<<","<<token.attr[2]<<") "<<token.attr[3]<<endl;
						cout<<"sidx = "<<token.attr[4]<<", Cid = "<<token.attr[5]<<endl;
						route.ShowPath(cb.GetwPath(token.attr[0]));
					}

					// analysis
					eval.AcceptEvent(token.attr[3],token.attr[4],route.Get_Path_Hop(cb.GetwPath(token.attr[0])));
				}
				else
				{
					if(show_detail)
					{
						cout<<"req "<<token.attr[0]<<" is blocked."<<endl;
						cout<<"("<<token.attr[1]<<","<<token.attr[2]<<") "<<token.attr[3]<<endl;
					}
					
					// analysis
					eval.BlockEvent(token.attr[3]);
				}
				
				// generate request
				Generate_request(G.node_num, ReqMaxSize, current_req_num, current_erlang);
				current_req_num++;

			} break;

			case FINISH:
			{
				if(show_detail)
					cout<<"req "<<token.attr[0]<<" FINISH"<<endl;

				// release op
				EOP_Release(&cb, token);
				update_completions();
			} break;

		}  // end of switch
	}  // end while

	// To process the remaining event in the schedule
	while(1)
	{
		int event_statu = next_event(&event,&token);  // get the event
		if(event_statu != ES_OK)
		{
			//cout<<"event_statu != ES_OK"<<endl;
			break;
		}
		switch(event)
		{
			case FINISH:
			{
				if(show_detail)
					cout<<"req "<<token.attr[0]<<" FINISH"<<endl;

				// release op
				EOP_Release(&cb, token);
				update_completions();
			} break;

			default:
			{
				//cout<<"default, event = "<<event<<endl;
			}
		}
	}

	// analysis
	eval.GetAnalsis(&ev);
	cb.GetAnalysis_1(&cba);
}

void Generate_request(int node_num, int ReqMaxSize, int req_id, int arrival_erlang)
{
	TOKEN token;

	// determine source, dest., FS size
	int source = rand() % node_num, dest;
	do
	{
		dest = rand() % node_num;
	} while(source == dest);

	// SLICE
	//int size = rand()%ReqMaxSize + 1;  // ­×¥¿
	int size = ceil(uniform(0.1,120)/8);  // 8 = capacity (2) * nod_level (4)

	/*
	// MLR
	int Mlen = 3;
	int MLR[3] = {4,10,24};
	double Mrate[3] = {100.0, 400.0, 1000.0};  // date rate
	int size = MLR[random(0,Mlen-1)];
	*/
	
	double arr = expntl(1.0/arrival_erlang);

	// set token
	token.attr[0] = req_id;  // the ID of the request.
	token.attr[1] = source;  // the source of req.
	token.attr[2] = dest;  // the dest. of req.
	token.attr[3] = size;  // the used FS size of req.
	token.attr[4] = token.attr[5] = -1;  // the starting index, cluster index. (may change in PROCESS)

	update_arrivals();
	schedule(PROCESS,arr,token);  // add event into schedule.
}

bool EOP_Setup(cbrsa* cb, TOKEN* token)
{
	// the request information can be got from the token
	int req_id = token->attr[0];
	int source = token->attr[1];
	int dest = token->attr[2];
	int size = token->attr[3];

	// counting time
	LARGE_INTEGER  large_interger;
	QueryPerformanceFrequency(&large_interger);
	double dff = large_interger.QuadPart;
	__int64 tp1, tp2;
	// start tp1
	QueryPerformanceCounter(&large_interger);
	tp1 = large_interger.QuadPart;

	//double tp1 = clock();
	// Finding a working cluster and a working path to assign request
	int ci = -1;  // the index of the cluster
	// run the corresponding method.
	switch(RS)
	{
		case KSP:
		{
			switch(CAS)
			{
				case SP:
				{
					ci = cb->FindCwToAssign_KSP_SP(req_id, size, route.Get_K_value(), route.ReturnKPath(source, dest));
					//ci = cb->FindCwToAssign_AR_SP_2(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));
				} break;
				case FF:
				{ 
					ci = cb->FindCwToAssign_KSP_FF(req_id, size, route.Get_K_value(), route.ReturnKPath(source, dest));
					//ci = cb->FindCwToAssign_AR_MF(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));
				} break;
				case MF:
				{
					//ci = cb->FindCwToAssign_KSP_MF(req_id, size, route.Get_K_value(), route.ReturnKPath(source, dest));
					ci = cb->FindCwToAssign_KSP_MF_1(req_id, size, route.Get_K_value(), route.ReturnKPath(source, dest));
					//ci = cb->FindCwToAssign_KSP_MF_2(req_id, size, route.Get_K_value(), route.ReturnKPath(source, dest));
				} break;
			}
		} break;
		case AR:
		{
			switch(CAS)
			{
				case SP:
				{
					//ci = cb->FindCwToAssign_AR_SP_2(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // SP
					//ci = cb->FindCwToAssign_AR_SP_3(req_id, size, source, dest, route.GetPathNum(source,dest), route.ReturnKPath(source, dest));
					//ci = cb->FindCwToAssign_AR_SP_4(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // SP_MF
					ci = cb->FindCwToAssign_AR_SP_6(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // SP_MF_ss
				} break;
				case FF:
				{
					ci = cb->FindCwToAssign_AR_FF(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));
					//ci = cb->FindCwToAssign_AR_MF_1(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));
					//ci = cb->FindCwToAssign_AR_SP_2(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));
				} break;
				case MF: 
				{
					//ci = cb->FindCwToAssign_AR_MF(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // MF
					//ci = cb->FindCwToAssign_AR_MF_1(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // MF_SP
					//ci = cb->FindCwToAssign_AR_MF_2(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // MC
					//ci = cb->FindCwToAssign_AR_MF_3(req_id, size, source, dest, route.GetPathNum(source,dest), route.ReturnKPath(source, dest));
					//ci = cb->FindCwToAssign_AR_MF_4(req_id, size, source, dest, route.GetPathNum(source,dest), route.ReturnKPath(source, dest));
					//ci = cb->FindCwToAssign_AR_MF_5(req_id, size, source, dest, route.GetPathNum(source,dest), route.ReturnKPath(source, dest));
					ci = cb->FindCwToAssign_AR_MF_6(req_id, size, source, dest, route.ReturnSinglePath(source, dest, 0));  // MC_SP

					//double add_time = route.Get_Path_Time(source, dest, cb->found_times+1);
					//rsa_time += add_time;
					//r_time += add_time;

				} break;
			}
		} break;
		default:
		{
			cout<<"No corresponding method."<<endl;
			exit(0);
		}
	}

	// end tp1
	QueryPerformanceCounter(&large_interger);
	tp1 = large_interger.QuadPart - tp1;
	r_time += tp1*1000/dff;

	// Defragmentation
	//if(ci==-1 && DF_flag)  // trigger the DF process
		//ci = cb->DefragProcess(req_id, size, RS, source, dest, route.Get_K_value(), route.ReturnKPath(source, dest));

	// if eop is acceptable.
	bool accept_eop = ci!=-1?true:false;
	if(accept_eop)  // if the EOP can be setup
	{
		// start tp2
		QueryPerformanceCounter(&large_interger);
		tp2 = large_interger.QuadPart;

		// assign req.
		//cb->AssignRequestIntoCw_1(req_id, size, ci);
		cb->AssignRequest(req_id);
		token->attr[4] = cb->Get_Cw_sidx(ci);
		token->attr[5] = ci;

		// end tp2
		QueryPerformanceCounter(&large_interger);
		tp2 = large_interger.QuadPart - tp2;
		sa_time += tp2*1000/dff;

		return true;
	}
	else  // else the EOP would be blocked
	{
		return false;
	}
}

void EOP_Release(cbrsa *cb, TOKEN token)
{
	int req_id = token.attr[0];  // the request's id
	//int source = token.attr[1];
	//int dest = token.attr[2];
	//int size = token.attr[3];
	//int sidx = token.attr[4];
	//int size = token.attr[3];
	//int widx = token.attr[5];  // the index of the working cluster

	//cb->Release_Request(req_id);
	//cb->Release_Request_1(req_id);
	cb->ReleaseRequest(req_id);
}

void Set_parameter(string* tf, int* K, int* Fmax, int* TotalReqNum, int* ReqMaxSize, int* IniErlang, int* IntervalErlang, int* ErlangTimes, int* EvalTimes)
{
	fstream file;
	string line;
	file.open(".\\input\\parameter.txt", ios::in);
	if(!file)
    {
        cout<<"No File parameter.txt"<<endl;
        exit(0);
    }
	else
	{
		while(!file.eof())
		{
			getline(file, line);
			int index = line.find(":");
			//cout<<line.substr(0,index)<<" = "<<line.substr(index+1,line.length()-1)<<endl;
			if(!line.substr(0,index).compare("Topology"))
				*tf = line.substr(index+1,line.length()-1);
			else if(!line.substr(0,index).compare("K"))
				istringstream(line.substr(index+1,line.length()-1)) >> *K;
			else if(!line.substr(0,index).compare("Fmax"))
				istringstream(line.substr(index+1,line.length()-1)) >> *Fmax;
			else if(!line.substr(0,index).compare("TotalReqNum"))
				istringstream(line.substr(index+1,line.length()-1)) >> *TotalReqNum;
			else if(!line.substr(0,index).compare("ReqMaxSize"))
				istringstream(line.substr(index+1,line.length()-1)) >> *ReqMaxSize;
			else if(!line.substr(0,index).compare("IniErlang"))
				istringstream(line.substr(index+1,line.length()-1)) >> *IniErlang;
			else if(!line.substr(0,index).compare("IntervalErlang"))
				istringstream(line.substr(index+1,line.length()-1)) >> *IntervalErlang;
			else if(!line.substr(0,index).compare("ErlangTimes"))
				istringstream(line.substr(index+1,line.length()-1)) >> *ErlangTimes;
			else if(!line.substr(0,index).compare("EvalTimes"))
				istringstream(line.substr(index+1,line.length()-1)) >> *EvalTimes;
			else
			{
				// nothing
			}
		}  // end while
	}  // end if
}

string GetMethodOuptputFile()
{
	string str = "CBRSA";
	//"CBRSA_KSP_SP.txt"
	//switch()
	switch(RS)
	{
		case KSP: str.append("_KSP"); break;
		case AR: str.append("_AR"); break;
	}
	switch(CAS)
	{
		case FF: str.append("_FF"); break;
		case SP: str.append("_SP"); break;
		case MF: str.append("_MF"); break;
		case Uvalue: str.append("_Uvalue"); break;
	}
	//str.append(".txt");
	return str;
}

void initEvalvalue(EvalValue* ev, int ReqMaxSize)
{
	ev->BBP = 0;
	ev->BP = 0;
	ev->MS = 0;
	ev->usedFS = 0;
	ev->duration = 0;
	ev->path_hop = 0;
	ev->bp_size.assign(ReqMaxSize+1, 0);
	ev->r_time = 0;
	ev->sa_time = 0;
}

void initCbAalysis(cb_Analysis* cba, int ReqMaxSize, int N)
{
	cba->r_count = 0;
	cba->size_block = 0;
	//----------
	cba->avg_cw_size = 0;
	cba->avg_continuous_cw_assign = 0;
}

Topology Get_topology(string file_name)
{
	Topology G;
    // Read Topology
    fstream file;
    char line[40];
    string filepath = ".\\input\\" + file_name;
    file.open(filepath.c_str(), ios::in);
    if(!file)
    {
        cout<<"No Topology File ./input/"<<file_name<<"exists."<<endl;
        exit(0);
    }
    else
    {
        // get the # of nodes
        file.getline(line, sizeof(line));
        istringstream(line) >> G.node_num;
		// tmp graph to store the info of topology.
		vector< vector<int> > g(G.node_num,vector<int>(G.node_num,0));
        // set value in "graph"
        char *split;
        int i = 0;
        while(!file.eof())
        {
            file.getline(line, sizeof(line));
            //cout << line << endl;
            split = strtok(line," ");
            int j = 0;
            while(split != NULL)
            {
                //cout << "aa" << split << endl;
                istringstream(split) >> g[i][j++];
                split = strtok(NULL," ");
            }
            i++;
        }
        file.close();
		// copy g into two dim array "graph"
		G.graph.assign(g.begin(), g.end());
    }
    //ShowTwoDimArray(G.node_num, G.node_num, G.graph);  // test
	return G;
}

void Set_method(int argc, char* argv[])
{
	if(argc == 1)  // default
	{
		RS = KSP;
		CAS = FF;
	}
	else
	{
		if(!strcmp(argv[1],"KSP"))
			RS = KSP;
		else if(!strcmp(argv[1],"AR"))
			RS = AR;

		if(!strcmp(argv[2],"FF"))
			CAS = FF;
		else if(!strcmp(argv[2],"SP"))
			CAS = SP;
		else if(!strcmp(argv[2],"MF"))
			CAS = MF;
	}
}
