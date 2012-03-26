#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

class Equation
{
public:
	int numParameters;
	vector<int> reactants;
	vector<int> products;
	vector<int> AIs;
	
	Equation(Equation* e); /**/
	Equation(vector<Equation*>, int index); /**/
	Equation(vector<int> reacts, vector<int> products); /**/
	~Equation();
	void generateActivatorsInhibitors(vector<bool> used); /**/
	string print();
	
	double getVValue(vector<double> concentrations, vector<double> K, double kCatMinus, double kCatPlus); /**/
}; 
