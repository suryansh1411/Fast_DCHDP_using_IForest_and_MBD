#include <iostream>
#include <queue>
#include "./data.cpp"
#include "distance.cpp"
#include "./iforest.cpp"
#include "./dchdp.cpp"


int find(vector<int>& par, int x)
{
    if(par[x]==x) return x;
    return par[x]=find(par, par[x]);
}

void myunion(vector<int>& par, int x, int y, vector<int>& rank)
{
    int par_x = find(par, x);
    int par_y = find(par, y);

    if(rank[par_x]==rank[par_y])
    {
        par[par_x]=par_y;
        rank[par_y]++;
    }
    else if(rank[par_x]>rank[par_y])
    {
        par[par_y]=par_x;
    }
    else
    {
        par[par_x]=par_y;
    }
}

vector<int> merge(vector<pair<int, int>> mergers, int numClasses)
{
    vector<int> par(mergers.size()+1);
    vector<int> rank(mergers.size()+1, 0);
    for(int i=0; i<par.size(); i++)
    {
        par[i]=i;
    }
    // cout<<mergers.size()<<endl;

    for(int i=0; i<mergers.size()+1-numClasses; i++)
    {
        myunion(par, mergers[i].first, mergers[i].second, rank);    
    }
    
    for(int i=0; i<par.size(); i++)
    {
        par[i]=find(par, i);
    }
    return par;
}




int main(int argc, char* argv[])      //(argv[1] = inputdataFile.csv , argv[2] = number of clusters
{
	srand(time(0));
//read input from csv file
	cout<<"Dc dPTime iforestTime find_markedNodesTime find_potential_dcNN_listTime relativedistanceTime ccidentificationTime_without_k ccidentificationTime clusterassignTime "<<endl;
	const string dataset = argv[1];
	
	struct timespec start_dP,end_dP;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_dP);

	data *dataObj = new data();
    dataObj->createDataVector(dataset+"/"+dataset+".csv");
    const data &refDataObject = *dataObj;
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_dP);
    double dPTime =  (((end_dP.tv_sec - start_dP.tv_sec) * 1e9)+(end_dP.tv_nsec - start_dP.tv_nsec))*1e-9;



//claculate iforest construction

	
	struct timespec start_iforest,end_iforest;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_iforest);
	vector<vector<double>> massMatrix;
	
	iforest *iForestObject = new iforest(refDataObject, 100, 256);
	iforest &refiForestObject = *iForestObject;
	
	iForestObject->constructiForest();
	iForestObject->computeNodeMass();
	iForestObject->computeLCA_lookup();
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_iforest);
    double iforestTime =  (((end_iforest.tv_sec - start_iforest.tv_sec) * 1e9)+(end_iforest.tv_nsec - start_iforest.tv_nsec))*1e-9;
	//cout<<"MatrixTime="<<iforestTime<<endl;
	
	//cout<<"iforest done"<<endl;
	iForestObject->write_smallest_largest_leaf(dataset);

//read number of cluster as a user input
	const double k = atof(argv[2]);
	//cout<<"k reading done"<<endl;



//clustering start 
	
	ofstream write_clusterCenters(dataset+"/intermediatefiles/smpc_clusterCenters.csv",ios::out|ios::binary);
	if(!write_clusterCenters){
		cout<<"Can not open input data file: intermediatefiles/smpc_clusterCenters.csv"<<endl;
		exit(0);
	}
	write_clusterCenters.close();

	ofstream write_cId(dataset+"/intermediatefiles/mpc_cId.csv",ios::out|ios::binary);
	if(!write_cId){
		cout<<"Can not open input data file: intermediatefiles/smpc_cId.csv"<<endl;
		exit(0);
	}
	write_cId<<"pointid actual predicted"<<endl;
    write_cId.close();



	//cout<<"both cid an dhaloid files initialized"<<endl;

	for(double Dc = 0.1; Dc < 1; Dc = Dc+0.1){
		int pointi = 0;
		DCHDP *dchdp = new DCHDP(refiForestObject, massMatrix, Dc);  

		struct timespec start_find_markedNodes,end_find_markedNodes;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_find_markedNodes);
		
	
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_find_markedNodes);
		double find_markedNodesTime =  (((end_find_markedNodes.tv_sec - start_find_markedNodes.tv_sec)* 1e9)+(end_find_markedNodes.tv_nsec - start_find_markedNodes.tv_nsec))*1e-9;

		struct timespec start_find_potential_dcNN_list,end_find_potential_dcNN_list;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_find_potential_dcNN_list);
		//density
		dchdp->find_potential_dcNN_list();

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_find_potential_dcNN_list);
		double find_potential_dcNN_listTime =  (((end_find_potential_dcNN_list.tv_sec - start_find_potential_dcNN_list.tv_sec) * 1e9)+(end_find_potential_dcNN_list.tv_nsec - start_find_potential_dcNN_list.tv_nsec))*1e-9;




	//claculate relative distance
		struct timespec start_relativedistance,end_relativedistance;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_relativedistance);
		
		//nndh
	
		dchdp->findDensityConnectivity();
		dchdp->sRelativeDistance();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_relativedistance);
    	double relativedistanceTime =  (((end_relativedistance.tv_sec - start_relativedistance.tv_sec) * 1e9)+(end_relativedistance.tv_nsec - start_relativedistance.tv_nsec))*1e-9;
	//cluster center identification with k

		
		struct timespec start_ccidentification_without_k,end_ccidentification_without_k;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_ccidentification_without_k);

		//mPeakObj->numberOfClusterCenters();

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_ccidentification_without_k);
    
		double ccidentificationTime_without_k =  (((end_ccidentification_without_k.tv_sec - start_ccidentification_without_k.tv_sec) * 1e9)+(end_ccidentification_without_k.tv_nsec - start_ccidentification_without_k.tv_nsec))*1e-9;

		struct timespec start_ccidentification,end_ccidentification;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_ccidentification);

		// mPeakObj->assignClusterId(k);

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_ccidentification);
    
		double ccidentificationTime =  (((end_ccidentification.tv_sec - start_ccidentification.tv_sec) * 1e9)+(end_ccidentification.tv_nsec - start_ccidentification.tv_nsec))*1e-9;
		//cout<<"cc identification done"<<endl;

	//cluster assignment
		struct timespec start_clusterassign,end_clusterassign;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_clusterassign);
	
		vector<pair<int, int>> mergers = dchdp->clusterAssignment();
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_clusterassign);
		double clusterassignTime =  (((end_clusterassign.tv_sec - start_clusterassign.tv_sec) * 1e9)+(end_clusterassign.tv_nsec - start_clusterassign.tv_nsec))*1e-9;
		
	// merge
		vector<int> clusters = merge(mergers, k);

		ofstream write_cId(dataset+"/intermediatefiles/mpc_cId.csv",ios::app|ios::binary);
		pointi = 0;
		for(auto i:clusters){
			write_cId<<pointi<<" "<<refDataObject.dataVector[pointi]->label<<" "<<i<<endl;
			pointi++;
		}
		write_cId.close();

		//cout<<"cid writing done"<<endl;
	
		//cout<<sizeof(massMatrix)<<"+"<<sizeof(mPeakObj)<<" ";
						
		delete dchdp;

		cout<<Dc<<" "<<dPTime<<" "<<iforestTime<<" "<<find_markedNodesTime<<" "<<find_potential_dcNN_listTime<<" "<<relativedistanceTime<<" "<<ccidentificationTime_without_k<<" "<<ccidentificationTime<<" "<<clusterassignTime<<endl;
	}































}
