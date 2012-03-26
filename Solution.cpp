//Solution.cpp and ODE Function
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "Solution.h"

string Solution::print()
{
	stringstream s;
	s << "Solution------------" << endl;
	s << "NumParams: " << numParameters << endl;
	for(int i=0; i < equations.size(); i++)
	{
		s << equations[i]->print();
	}
	s << "K values:" << endl;
	for(int i=0; i < K.size(); i++)
	{
		if (metaboliteUsed[i])
			s << K[i] << " ";
	}
	s << endl;
	s << "Kplus: " << kCatPlus << " Kmnius: " << kCatMinus << endl;
	if (errorValue != INT_MAX)
		s << "Error: " << errorValue << endl;
	return s.str();
}

string Solution::getKs()
{
	stringstream output;
	for(int i=0; i < metaboliteUsed.size(); i++)
	{
		if (metaboliteUsed[i])
		{
			output << K[i] << "\t";
		}
		else
		{
			output << "N\\A\t";
		}
	}
	output << kCatPlus << "\t" << kCatMinus << endl;
	return output.str();
}

bool Solution::operator<( const Solution& rhs )
{ 
	return (errorValue < rhs.errorValue); 
}

bool Solution::operator==(const Solution& rhs)
{
	for(int i=0; i < metaboliteUsed.size(); i++)
	{
		if (metaboliteUsed[i] != rhs.metaboliteUsed[i])
			return false;
	}
	return true;
}
	
Solution::Solution()/**/
{
	errorValue = INT_MAX;
}

Solution::Solution(int testCase)/**/
{
	errorValue = INT_MAX;
	switch (testCase)
	{
		case 1:
			numParameters = 6;
			parseInitialEquations(testCase);
			break;
		case 2:
			numParameters = 10;
			parseInitialEquations(testCase);
			break;
		case 3:
			numParameters = 19;
			parseInitialEquations(testCase);
	}
	for(int i=0; i < numParameters; i++)
	{
		if (testCase == 3) //preprocessing for larger networks (things we know will be combined)
		{
			if ( i == 0 || i == 1 || i == 2 || i == 3 || i == 12 || i == 13 || i == 16 || i == 18)
			{
				metaboliteUsed.push_back(true);
			}
			else
			{
				metaboliteUsed.push_back(false);
			}
		}
		else
		{	
			metaboliteUsed.push_back(true);
		}	
		K.push_back(0);
	}
	kCatPlus = 0;
	kCatMinus = 0;
	
}

Solution::Solution(Solution* s)/**/
{
	this->numParameters = s->numParameters;
	for(int i=0; i < s->equations.size(); i++)
	{
		this->equations.push_back(new Equation(s->equations[i]));
	}
	this->inputsOutputs = s->inputsOutputs;
	this->K = s->K;
	this->kCatPlus = s->kCatPlus;
	this->kCatMinus = s->kCatMinus;
	this->errorValue = s->errorValue;
	this->calculatedData = s->calculatedData;
	this->metaboliteUsed = s->metaboliteUsed;
	
}

Solution::~Solution()
{
	for(int i=0; i < equations.size(); i++)
	{
		delete equations[i];
	}
}

void Solution::parseInitialEquations(int testCase)/**/
{
	ifstream ifile;
	int x;
	switch (testCase)
	{
		case 1:
			ifile.open("Test1Params.txt");
			for(int j=0; j < 3; j++)
			{
				vector<int> reactants, products;
				char buffer[1000];
				ifile.getline(buffer, 1000);
				stringstream line(buffer);
				for(int i=0; i < numParameters; i++)
				{
					line >> x;
					reactants.push_back(x);
				}
				for(int i=0; i < numParameters; i++)
				{
					line >> x;
					products.push_back(x);
				}
				equations.push_back(new Equation(reactants, products));
			}
			inputsOutputs.push_back(true);
			inputsOutputs.push_back(false);
			inputsOutputs.push_back(false);
			inputsOutputs.push_back(true);
			inputsOutputs.push_back(true);
			inputsOutputs.push_back(true);
			
			ifile.close();
			
			break;
		case 2:
			ifile.open("Test2Params.txt");
			for(int j=0; j < 3; j++)
			{
				vector<int> reactants, products;
				char buffer[1000];
				ifile.getline(buffer, 1000);
				stringstream line(buffer);
				for(int i=0; i < numParameters; i++)
				{
					line >> x;
					reactants.push_back(x);
				}for(int i=0; i < numParameters; i++)
				{
					line >> x;
					products.push_back(x);
				}
				equations.push_back(new Equation(reactants, products));
			}
			for(int i=0; i < numParameters; i++)
			{
				if (i == 0 || i == 1)
				{
					inputsOutputs.push_back(true);
				}
				else
				{
					inputsOutputs.push_back(false);
				}
			}
			ifile.close();
			break;
		case 3:
			ifile.open("Test3Params.txt");
			for(int j=0; j < 5; j++)
			{
				vector<int> reactants, products;
				char buffer[1000];
				ifile.getline(buffer, 1000);
				stringstream line(buffer);
				for(int i=0; i < numParameters; i++)
				{
					line >> x;
					reactants.push_back(x);
				}for(int i=0; i < numParameters; i++)
				{
					line >> x;
					products.push_back(x);
				}
				equations.push_back(new Equation(reactants, products));
			}
			for(int i=0; i < numParameters; i++)
			{
				if (i == 0 || i == 2 || i == 3 || i == 18)
				{
					inputsOutputs.push_back(true);
				}
				else
				{
					inputsOutputs.push_back(false);
				}
			}
			ifile.close();
			break;
	}
}

int Solution::getNumMetabolitesUsed()
{
	int numMetabolitesUsed = 0;
	for(int i=0; i < metaboliteUsed.size(); i++)
	{
		if (metaboliteUsed[i] && !inputsOutputs[i])
		{
			numMetabolitesUsed++;
		}
	}
	return numMetabolitesUsed;
}
Solution* Solution::getSolution(int index)/**/
{
	Solution* newSoln = new Solution();
	newSoln->numParameters = this->numParameters;
	newSoln->metaboliteUsed = this->metaboliteUsed;
	newSoln->metaboliteUsed[index] = false;
	newSoln->inputsOutputs = this->inputsOutputs;
	newSoln->K = this->K;
	newSoln->kCatPlus = this->kCatPlus;
	newSoln->kCatMinus = this->kCatMinus;

	int i=0;
	vector<Equation*> combineEquations;
	for(i=0; i < this->equations.size(); i++)
	{
		if (this->equations[i]->reactants[index]  != 0 || this->equations[i]->products[index] != 0)
		{
			this->equations[i]->generateActivatorsInhibitors(metaboliteUsed);
			combineEquations.push_back(this->equations[i]);
		}
		else
		{
			newSoln->equations.push_back(this->equations[i]);
		}
	}
	
	if(combineEquations.size() != 0)
	{
		Equation* newEq = new Equation(combineEquations, index);
		newSoln->equations.push_back(newEq);
	}
	
	return newSoln;
	
}

vector<Solution*> Solution::generateNewSolutions()/**/
{
	vector<Solution*> newSolutions;
	for(int i=0; i < numParameters; i++)
	{
		Solution* temp;
		if (!inputsOutputs[i] && metaboliteUsed[i])
		{
			temp = getSolution(i);
			if (temp != NULL)
			{
				newSolutions.push_back(temp);
			}
		}
	}
	return newSolutions; 
}

void Solution::generateRandomK(double minK, double maxK)/**/
{
	K.clear();
	for(int i=0; i < numParameters; i++)
	{
		double x = (maxK-minK)*((rand()%1000+1)/1000.0f) + minK;
		K.push_back(x);
	}
	kCatPlus = (maxK-minK)*((rand()%1000+1)/1000.0f) + minK;
	kCatMinus = (maxK-minK)*((rand()%1000+1)/1000.0f) + minK;
}

void Solution::generateNewAIs()/**/
{
	for(int i=0; i < equations.size(); i++)
	{
		equations[i]->generateActivatorsInhibitors(metaboliteUsed);
	}
}

string Solution::getKey()
{
	stringstream s;
	for(int i = 0; i < metaboliteUsed.size(); i++)
	{
		if (metaboliteUsed[i])
		{
			s << i << ',';
		}
	}
	return s.str().c_str();
}


void Solution::getVs(vector<double> &x, vector<double> &dxdt)/**/
{
	for(int i=0; i < x.size(); i++)
	{
		dxdt[i] = 0;
	}
	for(int i=0; i < equations.size(); i++) //go through every equation
	{
		for(int j=0; j < equations[i]->numParameters; j++) //for each parameter in the equation
		{
			if (equations[i]->reactants[j] != 0) //if it is a reactant, there is a negative change
			{
				dxdt[j] -= equations[i]->getVValue(x, K, kCatMinus, kCatPlus);
			}
			if (equations[i]->products[j] != 0) //if it is a product, there is a positive change
			{
				dxdt[j] += equations[i]->getVValue(x, K, kCatMinus, kCatPlus);
			}
		}
	}
}
