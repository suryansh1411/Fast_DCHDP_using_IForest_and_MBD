#define DISTTYPE Euclidean
#ifndef DISTANCE_H
#define DISTANCE_H
//#include "data.h"
#pragma once
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>



class Distance
{
    public:
    
	Distance(const data &, const string &);                              
    Distance();
        
	virtual ~Distance();
	
	void computeDistanceMatrix(vector<vector<double>> & matrix);
	double euclideanDistance(int, int);

	void computeMassMatrix(vector<vector<double>> & matrix, int, float);


	private:
    
	const string _distType;
	const data & _dataObject;	
	
		
};


#endif 
