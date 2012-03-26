//Equation.cpp
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "Equation.h"
string Equation::print()
{
	stringstream s;
	s << "Printing Equation" << endl;
	s << "Reactants" << endl;
	for(int i=0; i < reactants.size(); i++)
	{
		s << reactants[i] << " ";
	}
	s << endl << "Products" << endl;
	for(int i=0; i < products.size(); i++)
	{
		s << products[i] << " ";
	}
	s << endl << "AIs" << endl;
	for(int i=0; i < AIs.size(); i++)
	{
		s << AIs[i] << " ";
	}
	s << endl;
	return s.str();
}

Equation::Equation(Equation* e) /**/
{
	this->numParameters = e->numParameters;
	this->reactants = e->reactants;
	this->products = e->products;
	this->AIs = e->AIs;
}

Equation::Equation(vector<int> reacts, vector<int> prods) /**/
{
	this->reactants = reacts;
	this->products = prods;
	this->numParameters = reacts.size();
}

Equation::Equation(vector<Equation*> inputs, int index) /**/
{
	
	this->reactants = inputs[0]->reactants;
	this->reactants[index] = 0;
	
	this->products = inputs[0]->products;
	this->products[index] = 0;
	
	this->AIs = inputs[0]->AIs;
	this->AIs[index] = 0;
	
	this->numParameters = reactants.size();
	
	
	for(int i=1; i < inputs.size(); i++)
	{
		for(int j=0; j < reactants.size(); j++)
		{
			if (j != index)
			{
				if (inputs[i]->reactants[j] != 0)
					reactants[j] += inputs[i]->reactants[j];
				if (inputs[i]->products[j] != 0)
					products[j] += inputs[i]->products[j];
				if (inputs[i]->AIs[j] != 0)
					AIs[j] += inputs[i]->AIs[j];
			}
		}
	}
	/*
	for(int i=0; i < reactants.size(); i++)
	{
		if ( reactants[i] != 0 && products[i] != 0) //collision fix
		{
			if (reactants[i] > products[i]) //item is specifically a reactant
			{
				reactants[i] = 1;
				products[i] = 0;
			}
			else //item is specifically a product, this will happen if they are the same value as well
			{
				products[i] = 1;
				reactants[i] = 0;
			}	
		}
	}*/
	
}

Equation::~Equation()
{
}

void Equation::generateActivatorsInhibitors(vector<bool> used) /**/
{
	AIs.clear();
	for(int i=0; i < reactants.size(); i++)
	{
		if (used[i])
		{
			int chance = rand()%100;
			if (chance < 20)
			{
				AIs.push_back(1);
			}
			else if (chance < 40)
			{
				AIs.push_back(-1);
			}
			else
			{
				AIs.push_back(0);
			}
		}
		else
		{
			AIs.push_back(0);
		}
	}
}

double Equation::getVValue(vector<double> concentrations, vector<double> K, double kCatMinus, double kCatPlus) /**/
{
	double V = 1;
	for(int i=0; i < AIs.size(); i++) //handle inhibitors and activators
	{
		if (AIs[i] < 0) //inhibitor
		{
			V *= (K[i]/(K[i] + concentrations[i]));
		}
		if (AIs[i] > 0) //activator
		{
			V *= (1 + concentrations[i]/K[i]);
		}
	}
	double reactantProduct=kCatPlus, reactantSumProduct = 1, productProduct = kCatMinus, productSumProduct = 1;
	for(int i=0; i < reactants.size(); i++)
	{
		if (reactants[i] != 0)
		{
			reactantProduct *= (concentrations[i]/K[i]);
			reactantSumProduct *= (1 + concentrations[i]/K[i]);
		}
		if (products[i] != 0)
		{
			productProduct *= (concentrations[i]/K[i]);
			productSumProduct *= (1 + concentrations[i]/K[i]);
		}
	}
	V *= (reactantProduct - productProduct)/(reactantSumProduct + productSumProduct -1);
	return V;
}
