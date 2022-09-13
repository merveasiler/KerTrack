// @author Merve Asiler

#include "CommonUtils.h"
#include "BaseMathOpUtils.h"
#include <iostream>

vector<string> split(string inputStr, string delimiter) {

	vector<string> pieces;
	string data = inputStr, subdata = "";

	while (data != "") {
		int ind = data.find_first_of(delimiter);
		
		if (ind < 0)
			subdata = data;
		else
			subdata = data.substr(0, ind);
		
		if (subdata != "")
			pieces.push_back(subdata);

		if (ind < 0)
			break;
		data = data.substr(ind + 1);
	}

	return pieces;

}

tuple<vector<vector<int>>, vector<vector<double>>> kMeans(vector<double> numberSet, int k) {

	vector<vector<double>> clusterVals;
	vector<vector<int>> clusterKeys;
	vector<double> centers;
	clusterVals.resize(k);
	clusterKeys.resize(k);
	
	// randomly initialize centers
	for (int i = 0; i < k; i++) {
		int ind = rand() % numberSet.size();
		if (find(centers.begin(), centers.end(), numberSet[ind]) == centers.end())
			centers.push_back(numberSet[ind]);
		else
			i--;
	}

	bool is_repeat = true;
	int counter = 0;

	while (is_repeat) {
		for (int i = 0; i < k; i++) {
			clusterVals[i].clear();
			clusterKeys[i].clear();
		}

		for (int i = 0; i < numberSet.size(); i++) {
			double minDistance = numeric_limits<double>::infinity();
			int clusterInd;
			for (int j = 0; j < k; j++) {
				if (abs(numberSet[i] - centers[j]) < minDistance) {
					minDistance = abs(numberSet[i] - centers[j]);
					clusterInd = j;
				}
			}
			clusterVals[clusterInd].push_back(numberSet[i]);
			clusterKeys[clusterInd].push_back(i);
		}

		is_repeat = false;
		counter++;
		cout << counter << endl;

		for (int i = 0; i < k; i++) {
			double sum = 0;
			for (int j = 0; j < clusterVals[i].size(); j++)
				sum += clusterVals[i][j];
			double newCenter = sum / clusterVals[i].size();
			if (abs(newCenter - centers[i]) >= EPSILON)
				is_repeat = true;
			centers[i] = newCenter;
		}
	}

	tuple< vector<vector<int>>, vector<vector<double>> > clusters = make_tuple(clusterKeys, clusterVals);
	return clusters;
}


