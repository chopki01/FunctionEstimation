//Driver.cpp

#include <iostream>
#include "NetworkSimplifier.h"

using namespace std;

int main(int argc, char* argv[])
{
	int exNum;
	if (argc != 2)
	{
		cout << "Usage: " << argv[0] << " Example #\n";
		exit(1);
	}
	
	exNum = atoi(argv[1]);
	
	if (exNum == 0)
	{
		cout << "Invalid Example Number" << endl;
		exit(2);
	}
	
	NetworkSimplifier* simplifier = new NetworkSimplifier(exNum);
	cout << "Network Simplifier Created" << endl;
	simplifier->run();
	cout << "Network Simplifier Run" << endl;
	simplifier->outputResults();
	cout << "Finished Outputting Results" << endl;
	delete simplifier;
	exit(0);
}
