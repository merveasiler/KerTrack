// @author Merve Asiler

#pragma once

#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

vector<string> split(string inputStr, string delimiter);

tuple<vector<vector<int>>, vector<vector<double>>> kMeans(vector<double> numberSet, int k);