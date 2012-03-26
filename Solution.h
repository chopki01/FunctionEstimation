#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "Equation.h"

class Solution
{
  private:
	int numParameters;
	void parseInitialEquations(int index); /**/
	
  public:
    bool operator<( const Solution& rhs );
	bool operator==(const Solution& rhs );
	
	vector<Equation*> equations;
	vector<bool> metaboliteUsed;
	vector<bool> inputsOutputs;
	vector<double> K;
	double kCatPlus, kCatMinus;
	
	double errorValue; 
	vector<vector<vector<double> > > calculatedData; 

	Solution(); 			/**/
	Solution(int testCase); /**/
	Solution(Solution* s);  /**/
	~Solution(); 			/**/
	
	vector<Solution*> generateNewSolutions(); /**/
	void getVs(vector<double> &x, vector<double> &dxdt); /**/
	void generateRandomK(double minK, double maxK); /**/
	void generateNewAIs(); /**/
	string getKs(); //Returns string with k values
	
	int getNumMetabolitesUsed();
	string getKey();
	Solution* getSolution(int index);
	
	string print();
	
	
};
