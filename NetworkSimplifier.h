//SA.h

#include "Solution.h"
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <map>

class mycomparison
{
    bool reverse;
	public:
    mycomparison(){}
    bool operator() (Solution* lhs, Solution* rhs) 
    {
		if (lhs->errorValue != lhs->errorValue) //if lhs is nan, always return false
		{
			return false;
		}
		else
		{
			if (rhs->errorValue != rhs->errorValue) //if rhs nan and lhs not, return true
			{
				return true;
			}
			else //do actual comparison
			{
				return (*lhs < *rhs);
			}
		}
    }
};

class NetworkSimplifier
{
private:
 
	int exampleNumber;
	int numVars, numInputsOutputs;
	double maxK, minK;
	vector<bool> inputsOutputs;
	vector<string> metaboliteNames;
	
	
	vector< vector<double> > maxValues;
	vector< vector<double> > minValues;
	
  
	
	Solution* bestSolution;
	vector<Solution*> bestSolutions;
	vector<vector<Solution*> > allSolutions;
	map<string,Solution*> attemptedSolutions;
	int populationSize;
	vector<Solution*> population;
	
	vector<vector<vector<double> > > measuredData;
	vector<vector<vector<double> > > calculatedData;
	
	
	vector<vector<double> > parseFile(string filename); /**/
	void convertToCSV(string filename, vector<vector<double> > data);/**/
	void parseData();/**/
	void parseNames();
	void parseTestingData();
	void setMaxMinValues();
	void generateFigure1(vector<Solution*> bestSolns, vector<Solution*> worstSolns);
	void generateFigure2(vector<vector<double> > d);
	void generateFigure3();
	
	
	
	Solution* crossover(int index1, int index2);
	Solution* generationalCrossover(Solution* s, Solution* t);
	double evaluate();/**/
	bool callODESolver(); /**/
	void generateNewInitialPopulation(Solution* s);
	
	double calculateError();
	
	
	
public:
	
	
	NetworkSimplifier(int exampleNumber); /**/
	~NetworkSimplifier();
	void run(); /**/
	void outputResults();
	void outputTestingResults(vector<Solution*>);
};
