//SA.cpp
#include <climits>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "boost/numeric/odeint.hpp"
#include "NetworkSimplifier.h"

#define ITERATIONNUMBER 200
using namespace boost::numeric::odeint;
Solution* currentSolution;


void odes(vector< double > &x , vector< double > &dxdt , double t){
   currentSolution->getVs(x,dxdt);
}

NetworkSimplifier::NetworkSimplifier(int ex)/**/
{
	srand(time(NULL));
	populationSize = 500;
	this->exampleNumber = ex;
	Solution* firstSoln = new Solution(ex);
	numInputsOutputs = 0;
	
	
	switch (ex)
	{
		case 1:
			numVars = 6;
			maxK = 5;
			minK = 0;
			numInputsOutputs = 4;
			break;
		case 2:
			numVars = 10;
			maxK = 5;
			minK = 0;
			numInputsOutputs = 2;
			break;
		case 3:
			numVars = 19;
			maxK = 5;
			minK = 0;
			numInputsOutputs = 4;
			break;
	}
	
	for(int i=0; i< numVars; i++) //initiate starting array
	{
		vector<Solution* > solns;
		allSolutions.push_back(solns);
	}
	cout << "Metabolites Used: " << firstSoln->getNumMetabolitesUsed() << endl;
	allSolutions[firstSoln->getNumMetabolitesUsed()].push_back(firstSoln);
	
	for(int i=0; i < allSolutions.size(); i++)
	{
		cout << allSolutions[i].size() << " ";
	}
	cout << endl;
	
	inputsOutputs = firstSoln->inputsOutputs;
	attemptedSolutions.insert(pair<string,Solution*>(firstSoln->getKey(), firstSoln));
	
	parseData();
	parseNames();
	setMaxMinValues();
	
}

NetworkSimplifier::~NetworkSimplifier()
{
	if (bestSolution != NULL)
		delete bestSolution;
	for(int i=0; i < bestSolutions.size(); i++)
	{
		delete bestSolutions[i];
	}
	
}

void NetworkSimplifier::setMaxMinValues()
{
	for(int k=0;  k < measuredData.size(); k++)
	{
		vector<double> minInternal;
		vector<double> maxInternal;
		maxValues.push_back(maxInternal);
		minValues.push_back(minInternal);
		
		for(int i=0; i < measuredData[k][0].size(); i++) //isolate variable
		{
		  double max = INT_MIN, min = INT_MAX; //generate max and min values for each data set
		  for(int j=0; j < measuredData[k].size(); j++)
			{
			  max = (max < measuredData[k][j][i]) ? measuredData[k][j][i] : max;
			  min = (min > measuredData[k][j][i]) ? measuredData[k][j][i] : min;
			}
			maxValues[k].push_back(max);
			minValues[k].push_back(min);
		}
	}
}

void NetworkSimplifier::run() /**/
{
	int index = allSolutions.size()-1;
	while (index >= 0) //Run for each combination of solutions, removing each one etc
	{
		for(int i=0; i < allSolutions.size(); i++)
		{
			cout << allSolutions[i].size() << " ";
		}
		cout << endl;
	
		if (allSolutions[index].size() == 0 )
		{
			index--;
		}
		else
		{
			for(int i=0; i < allSolutions[index].size(); i++)
			{
				cout << "Generting New Initial Population" << endl;
				generateNewInitialPopulation(allSolutions[index][i]);
				cout << "Finished New Population" << endl;
				
				int iterations = 0;
				while (iterations++ < ITERATIONNUMBER /*Termination Conditions here */ ) //this is the amount of times to run GA on the current Solution
				{
					if (iterations % (ITERATIONNUMBER/4) == 0)
					{
						cout << "Currently on Iteration " << iterations << endl;
					}
					for(int j=0; j < populationSize; j++)
					{
						currentSolution = population[j];
						currentSolution->errorValue = evaluate();
						currentSolution->calculatedData = this->calculatedData;
					}
					
					sort(population.begin(), population.end(), mycomparison());	
					
					
					
					if (bestSolution == NULL || bestSolution->errorValue != bestSolution->errorValue)
					{	
						bestSolution = new Solution(population[0]); //make sure to create new solution so that there are no problems with deletions
						cout << "Best Solution Created at iteration " << iterations <<  endl;
						
						//bestSolution->print();
					}
					else if (population[0]->errorValue < bestSolution->errorValue && population[0]->errorValue == population[0]->errorValue)
					{
						cout << "New Best Solution at iteration " << iterations << endl;
						delete bestSolution;
						bestSolution = new Solution(population[0]); //ensures you don't delete memory twice
					}
					
					if (iterations == ITERATIONNUMBER) //if this was our last iteration break
						break;
						
					vector<Solution*> newPop;
					for(int j=0; j< population.size()/5; j++) //take top 20% of old population and keep them ##Elitism
					{
						newPop.push_back(new Solution(population[j]));
					}
					
					int choiceRange = population.size() * 4 / 5; //drop items 800-1000 ##Elitism
					
					while (newPop.size() != population.size()) //crossover/mutate rest of population
					{
						if (rand()%100 < 40 ) //20% chance for mutation
						{
							Solution* s = population[rand()%choiceRange]; //randomly choose one of the solutions
							s->generateNewAIs(); 			//mutate the solution
							s->generateRandomK(minK, maxK); //mutate solution
							newPop.push_back(new Solution(s));
						}
						else //80% chance for breeding
						{
							Solution* s = crossover(rand()%choiceRange, rand()%choiceRange);
							newPop.push_back(s);
						}
					}
					for(int i=0; i< population.size(); i++)//delete old memory
					{
						delete population[i];
					}
					population.clear();
					population = newPop;
				}
				cout << "Min Error Found: " << bestSolution->errorValue << endl;
				bestSolution->print();
				allSolutions[index][i] = bestSolution;
				
				//CLEANUP
				bestSolution = NULL;
				for(int j=0; j< population.size(); j++)//delete old memory
				{
					delete population[j];
				}
				population.clear();
			}
			
			//PostCondition: All Solutions in AllSolutions[index] are the best solutions for their specific number of params used
			//We now order them by their score, and generate the next generation from that 
			sort(allSolutions[index].begin(),allSolutions[index].end(), mycomparison());
			
			//generate breeding coefficients
			for(int i = 0; i < allSolutions[index].size(); i++) //only breed from top 80%
			{
				cout << "Now breeding for next generation. On index: " << i << endl;
				if (i <= 1) //elitism criteria
				{
					vector<Solution*> newSolns = allSolutions[index][i]->generateNewSolutions();
					cout << "Generated " << newSolns.size() << " new Solutions.\n";
					
					for(int j=0; j < newSolns.size(); j++)
					{
						if (attemptedSolutions.insert(pair<string,Solution*>(newSolns[j]->getKey(),newSolns[j])).second) //if this value was successfully inserted
						{
							allSolutions[newSolns[j]->getNumMetabolitesUsed()].push_back(newSolns[j]);
						}
						else
						{
							cout << "Found Duplicate" << endl;
						}
					}
				}
				else
				{
					//breed specific solutions
					float rate = 1.0f - ((float)i)/((float)allSolutions[index].size()); 
					for(int j=0; j < allSolutions[index].size(); j++)
					{
						if ( (rand()%100) < 20) //EDIT: Control crossover rate here
						{
							Solution* x = generationalCrossover(allSolutions[index][i], allSolutions[index][j]);
							if (attemptedSolutions.insert(pair<string,Solution*>(x->getKey(),x)).second) //if this value was successfully inserted
							{
								allSolutions[x->getNumMetabolitesUsed()].push_back(x);
							}
							else
							{
								cout << "Found Duplicate" << endl;
							}
						}
					}
				}
			}
			
			index--;
		}
		//END TODO
		
	}
	cout << "Finished" << endl;
}

Solution* NetworkSimplifier::crossover(int index1, int index2)
{
	Solution* s = new Solution(population[index1]);
	Solution* s2 = population[index2];
	
	
	for(int j=0; j < s->equations.size(); j++)
	{
		for(int i=0; i < s->equations[j]->AIs.size(); i++)
		{
			if (rand()%100 > 50)
			{
				s->equations[j]->AIs[i] = s2->equations[j]->AIs[i];
			}
		}
	}	
		
	for(int i=0; i < s->K.size(); i++)
	{
		if (rand()%100 > 50)
			s->K[i] = s2->K[i];
	}
	
	if (rand()%100 > 50)
		s->kCatPlus = s2->kCatPlus;
	if (rand()%100 > 50)
		s->kCatMinus = s2->kCatMinus;
		
	return s;
}

Solution* NetworkSimplifier::generationalCrossover(Solution* s, Solution* t)
{
	Solution* newSoln = new Solution(s);
	for(int i=0; i < s->metaboliteUsed.size(); i++)
	{
		if (!t->metaboliteUsed[i])
		{
			Solution* x = newSoln->getSolution(i); //get the solution by removing metabolite i
			//delete newSoln;
			newSoln = x;
		}
	}
	return newSoln;
}

double NetworkSimplifier::calculateError()
{
    double error = 0;
	for(int k=0; k < measuredData.size(); k++) //per data file
	{
		for(int i=0; i < measuredData[k][0].size(); i++) //number of snapshots total
		{
		  double internalError = 0; //calculated normalized error per varaible
			if (inputsOutputs[i])
			{
				for(int j=0; j < measuredData[k].size(); j++)
				{
					internalError += pow((measuredData[k][j][i] - calculatedData[k][j][i]), 2);
				}
				internalError = pow(internalError, 0.5);
				if (maxValues[k][i] != maxValues[k][i])
				{
					internalError = internalError/(maxValues[k][i]-minValues[k][i]);
				}
				else
				{
					internalError = internalError/maxValues[k][i];
				}
				
				error += internalError; //add normalized error per variable to total error
			}
		}
	}
    return (double)error/(measuredData.size()*measuredData[0].size()*numInputsOutputs); //properly scale the data
}

//Evaluate current Solution
double NetworkSimplifier::evaluate()/**/
{
    if (callODESolver())
    {
		return calculateError();
    }
    else
	{
        return NAN;
	}
}

bool NetworkSimplifier::callODESolver()/**/
{
	double t = 0.0;
	double dt = 0.1;
	double tf = 10;
	int i, j, s;
	bool nonNan = true;
	for(int k=0; k < calculatedData.size(); k++)
	{
		bool internalNonNan = true;
		stepper_euler< vector<double> > stepper;
		stepper.adjust_size( calculatedData[k][0] );
		for( t=0, i=0 ; t <= tf ; t += dt, i++ )
		{
		    for(j = 0; j < calculatedData[k][i].size(); j++)
			{
			   calculatedData[k][i+1][j] = calculatedData[k][i][j];
			}
			stepper.do_step( odes , calculatedData[k][i+1] , t , dt );
			//NAN Checking
			for(j = 0; j < calculatedData[k][i].size(); j++)
			{
				if (calculatedData[k][i][j] != calculatedData[k][i][j]) //shortcuts out if solving the system of equations returns NaN as a value
				{
					nonNan = false;
					internalNonNan = false;
					break;
				}
			}
			if (!internalNonNan)
				break;
		}
	}
	return nonNan;
}

void NetworkSimplifier::parseData()/**/ //parses in all data for experiment
{
	stringstream paramFile;
	paramFile << "Test" << exampleNumber << "Params.txt";
	for(int i = 1; i < 6; i++)
	{
		stringstream name;
		name << "Ex" << exampleNumber << "Set" << i << ".dat";
		
		vector<vector<double> > temp = parseFile(name.str());
		
		measuredData.push_back(temp);
		
		convertToCSV(name.str(), temp);
	}
	calculatedData = measuredData;
}
void NetworkSimplifier::parseNames()
{
	ifstream nameFile;
	stringstream name;
	
	name << "Ex" << exampleNumber << "DataLabels.txt";
	
	nameFile.open(name.str().c_str());
	for(int i=0; i < numVars; i++)
	{
		string temp;
		nameFile >> temp;
		metaboliteNames.push_back(temp);
	}
	
	nameFile.close();
}
void NetworkSimplifier::parseTestingData()
{
	stringstream paramFile;
	vector<vector<vector<double> > > x;
	measuredData =x;
	paramFile << "Test" << exampleNumber << "Params.txt";
	for(int i = 6; i < 11; i++)
	{
		stringstream name;
		name << "Ex" << exampleNumber << "Set" << i << ".dat";
		
		vector<vector<double> > temp = parseFile(name.str());
		
		measuredData.push_back(temp);
		
		convertToCSV(name.str(), temp);
	}
	calculatedData = measuredData;
}
vector< vector<double> > NetworkSimplifier::parseFile(string filename)/**/
{
    ifstream input;
    input.open(filename.c_str());
    
    vector< vector<double> > measuredData;
    char buffer[1000];
    int lines = 102;
    while(lines--)
    {
        input.getline(buffer, 1000);
        stringstream line(buffer);
        double time;
        line >> time;
        vector<double> snapshot;
        
        for(int i=0; i < numVars; i++)
        {
            double temp;
            line >> temp;
            snapshot.push_back(temp);
        }
        measuredData.push_back(snapshot);
        
    }
    return measuredData;
}

void NetworkSimplifier::convertToCSV(string filename, vector<vector< double> > x)/**/ //converts avector to a csv file
{
	ofstream output;
	output.open((filename + ".csv").c_str());
	for(int i=0; i <x.size(); i++)
    {
		stringstream line;
		for(int j=0; j < x[i].size(); j++)
			{
				if (j != 0)
					line << ", ";
				line << x[i][j];
			}
		output << line.str() << endl;
    }
	output.close();
}


void NetworkSimplifier::generateNewInitialPopulation(Solution* s)/**/
{
	population.clear();
	for(int i=0; i < populationSize; i++)
	{
		population.push_back(new Solution(s));
		population[i]->generateRandomK(minK,maxK);
		population[i]->generateNewAIs();
	}
}

void NetworkSimplifier::outputTestingResults(vector<Solution*> bestSolutions)
{
  cout << "Parsing Testing Results\n";
  parseTestingData();
  cout << "Outputting Testing Results\n";


  for(int i=0; i< bestSolutions.size(); i++)
    {
      if (bestSolutions[i] != NULL)
		{
		  currentSolution = bestSolutions[i];
		  currentSolution->errorValue = evaluate();
		  cout << "Error for " << i <<  " Internal metabolites: " << currentSolution->errorValue << endl;
		  currentSolution->calculatedData = calculatedData;


		  for(int m=1; m < 6; m++)
			{
			  stringstream name;
			  ofstream finalData;
			  name << "Ex" << exampleNumber << "Set" << m+5 << "ResultsFor" << bestSolutions[i]->getNumMetabolitesUsed()  << "MetabolitesUsed.csv";
			  finalData.open(name.str().c_str());
			  cout << "File Opened" << endl;

			  for(int j=0; j < calculatedData[m-1].size(); j++)
				{
				  stringstream line;
				  for(int k=0; k < calculatedData[m-1][j].size(); k++)
					{
					  if (bestSolutions[i]->inputsOutputs[k])
					{
					  if (k != 0){ line << ", "; }
					  line << calculatedData[m-1][j][k];
					  line << ", ";
					  line << measuredData[m-1][j][k];
					}
					}
				  finalData << line.str() << endl;
				}
			  finalData.close();
			}
		}
    }
}
void NetworkSimplifier::generateFigure1(vector<Solution*> bestSolns, vector<Solution*> worstSolns)
{
	stringstream name;
	name << "Ex" << exampleNumber << "Figure1.csv";
	ofstream fig1;
	fig1.open(name.str().c_str());
	for(int i=0; i < bestSolns.size(); i++)
	{
		fig1 << bestSolns[i]->errorValue << ", " << worstSolns[i]->errorValue << endl;
	}
	fig1.close();
}
void NetworkSimplifier::generateFigure2(vector<vector<double> > d)
{
	stringstream name;
	name << "Ex" << exampleNumber << "Figure2.csv";
	ofstream fig2;
	fig2.open(name.str().c_str());
	for(int i=0; i < d.size(); i++)
	{
		sort(d[i].begin(), d[i].end());
	}
	for(int i=0; i < d.size(); i++)
	{
		for(int j=0; j < d[i].size(); j++)
		{
			if (j != 0)
				fig2 << ", ";
			fig2 << d[i][j];
			
		}
		fig2 << endl;
	}
	fig2.close();
}
void NetworkSimplifier::generateFigure3()
{
	ofstream fig3;
	stringstream name;
	name << "Ex" << exampleNumber << "Figure3.csv";
	fig3.open(name.str().c_str());
	for(int i=0; i < metaboliteNames.size(); i++)
	{
		fig3 << "k" << metaboliteNames[i] << "\t";
	}
	fig3 << "kCatPlus\tkCatMinus\n";
	
	for(int i=0; i < allSolutions.size(); i++)
	{
		for(int j=0; i < allSolutions[i].size(); i++)
		{
			fig3 << allSolutions[i][j]->getKs();
		}
	}
	
	fig3.close();
}

void NetworkSimplifier::outputResults()
{
	cout << "Outputting Results" << endl;
	Solution* finalSolution;
	double finalError = INT_MAX;
	vector<Solution*> bestSolutions;
	vector<Solution*> worstSolutions;
	
	
	vector<vector<double> > errorPerMetabolite; //Setup for figure 2
	for(int i=0; i < allSolutions.size(); i++)
	{
		vector<double> d;
		errorPerMetabolite.push_back(d);
	}
	
	
	ofstream params;
	ofstream eqns;
	eqns.open("eqns.txt");
	params.open("Params.txt");
	for(int i=0; i < allSolutions.size(); i++)
	{
		Solution* bestSoln = NULL;
		Solution* worstSoln = NULL;
		double error = INT_MAX;
		double worstError =-1;
		for(int j=0; j < allSolutions[i].size(); j++)
		{
			if (allSolutions[i][j]->errorValue < error)
			{
				bestSoln = allSolutions[i][j];
				error = bestSoln->errorValue;
			}
			if (allSolutions[i][j]->errorValue > worstError)
			{
				worstSoln = allSolutions[i][j];
				worstError = worstSoln->errorValue;
			}
			
			if (allSolutions[i][j]->errorValue < finalError)
			{
				finalSolution = allSolutions[i][j];
				finalError = finalSolution->errorValue;
			}
			
			for(int k=0; k < allSolutions[i][j]->metaboliteUsed.size(); k++)
			{
				if (allSolutions[i][j]->metaboliteUsed[k])
					errorPerMetabolite[k].push_back(allSolutions[i][j]->errorValue);
			}
		}
		if (bestSoln != NULL)
			bestSolutions.push_back(bestSoln);
		if (worstSoln != NULL)
			worstSolutions.push_back(worstSoln);
		
	}
	generateFigure1(bestSolutions,worstSolutions);
	generateFigure2(errorPerMetabolite);
	for(int i=0; i < bestSolutions.size(); i++)
	{
		if (bestSolutions[i] != NULL)
		{
			eqns << "______________________________________\n";
			eqns << "Best Solution for " << bestSolutions[i]->getNumMetabolitesUsed() <<  " internal metabolites:" <<  endl;
			eqns << bestSolutions[i]->print();
			eqns << endl << endl << "______________________________________\n";
			
			calculatedData = bestSolutions[i]->calculatedData;
			
			params << i << " internals used" << endl;
			for(int m=0; m < bestSolutions[i]->K.size(); m++)
			{
				params << bestSolutions[i]->K[m] << " ";
			}
			params << bestSolutions[i]->kCatPlus << " " << bestSolutions[i]->kCatMinus << endl;
			
			
			for(int m=1; m < 6; m++)
			{
				stringstream name;
				ofstream finalData;
				name << "Ex" << exampleNumber << "Set" << m << "ResultsFor" << bestSolutions[i]->getNumMetabolitesUsed()  << "MetabolitesUsed.csv";
				finalData.open(name.str().c_str());
				cout << "File Opened" << endl;
				
				for(int j=0; j < calculatedData[m-1].size(); j++)
				{
					stringstream line;
					for(int k=0; k < calculatedData[m-1][j].size(); k++)
					{
						if (bestSolutions[i]->inputsOutputs[k])
						{
							if (k != 0){ line << ", "; }
							line << calculatedData[m-1][j][k];
							line << ", ";
							line << measuredData[m-1][j][k];
						}
					}
					finalData << line.str() << endl;
				}
				finalData.close();
			}
		}
	}
	params.close();
	eqns.close();
	outputTestingResults(bestSolutions);
	
}
