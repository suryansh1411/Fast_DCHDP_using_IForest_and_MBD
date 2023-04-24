#ifndef mPeak_H
#define mPeak_H
#pragma once
#include "data.h"
#include "distance.h"
#include "data.cpp"
#include "distance.cpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>

class DCHDP
{
    public:
        
    DCHDP( iforest &,vector<vector<double>> &, const double &);
	//mPeakClustering( iforest &, const double &);
	virtual ~DCHDP();
	void find_markedNodes(treenode *, vector<vector<char>> &, int , int );
	void find_potential_dcNN_list();
	void computeLocalMass();
	void computeRelativeDistance();
	void sRelativeDistance();
	int numberOfClusterCenters();
	void assignClusterId(int k);
	void developCluster(int , int );
	void find_k_clustercenters(int);
	vector<pair<int, int>> clusterAssignment();
	void findDensityConnectivity();
	vector<double> getDeltaVector() const;
    vector<int> getNeighbourVector() const;
    vector<pair<double,int>>getOrdRho();
   // private: 
    
	vector<vector<double>> & _massMatrix;
	iforest & iForestObject;
	double _Dc;
	int _totalPoints;
	double _maxDistance;
	//vector<unordered_map<int,double>> _potential_dcNN_list;
    vector<vector<pair<double,int>>> _potential_dcNN_list;
    vector<double> _delta;
	vector<double> _rho;
	vector<pair<double,int>>_rhodelta;
    vector<int> _nneigh;
    vector<pair<double,int>> _ordrho;
	vector<int> _clusterCenters;
	vector<vector<int>> _clusters;
	vector<int> _cId;
	vector<int> _density_connectivity;
};

#endif 
