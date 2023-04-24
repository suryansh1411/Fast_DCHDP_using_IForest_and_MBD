#define DISTTYPE Euclidean
#include "distance.h"
#pragma once
#include "data.h"
#include "./iforest.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>

Distance::~Distance(){}
Distance::Distance(const data & dataObject, const string & distType): _dataObject(dataObject),_distType(distType){}

// method to compute distance between two points based on the provided distance type 
double Distance::euclideanDistance(int point1, int point2){
    double ans = 0;
		int dimensions = _dataObject.getnumAttributes();
        for(int i = 0; i<dimensions;i++){
            ans+= pow((_dataObject.dataVector[point1]->attributes[i] - _dataObject.dataVector[point2]->attributes[i]),2);
        }
        ans = sqrt(ans);
    return ans;
}


void Distance::computeDistanceMatrix(vector<vector<double>> & matrix){
    int size = _dataObject.getnumInstances();
    matrix.resize(size);
    for(int pointi = 0; pointi < size; pointi++){
        for(int pointj = 0; pointj < pointi; pointj++){
            matrix[pointi].push_back(euclideanDistance(pointi, pointj));
        }
    }
}

void Distance::computeMassMatrix(vector<vector<double>> & matrix, int numiTrees, float samplingFactor){
	iforest *iForestObject = new iforest(_dataObject, numiTrees, samplingFactor);
	iForestObject->constructiForest();
	iForestObject->computeNodeMass();
	iForestObject->computeMassMatrix(matrix); 
	delete iForestObject;
}

























