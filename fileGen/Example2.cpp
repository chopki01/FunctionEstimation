#include <iostream>
#include <iterator>
#include <boost/numeric/odeint.hpp>
#include <boost/lambda/lambda.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

typedef vector< double > concentrations;

#define PARAMETER_COUNT 36
#define VARIABLE_COUNT 10
double parameter [PARAMETER_COUNT];

void initParameters(double* param){
   param[0] = 1;
   param[1] = 1;
   param[2] = 2;
   param[3] = 1;
   param[4] = 2;
   param[5] = 1;
   param[6] = 1;
   param[7] = 1;
   param[8] = 2;
   param[9] = 1;
   param[10] = 2;
   param[11] = 1;
   param[12] = 1;
   param[13] = 1;
   param[14] = 2;
   param[15] = 1;
   param[16] = 2;
   param[17] = 1;
   param[18] = 0.1;
   param[19] = 1;
   param[20] = 0.1;
   param[21] = 0.1;
   param[22] = 1;
   param[23] = 0.1;
   param[24] = 0.1;
   param[25] = 1;
   param[26] = 0.1;
   param[27] = 1;
   param[28] = 1;
   param[29] = 1;
   param[30] = 1;
   param[31] = 1;
   param[32] = 1;
   param[33] = 1;
   param[34] = 1;
   param[35] = 1;
}

void initVariables(concentrations &variables){
	variables[0] = .5;
	variables[1] = .5;
	variables[2] = 1;
	variables[3] = 1;
	variables[4] = 1;
	variables[5] = 1;
	variables[6] = 1;
	variables[7] = 1;
	variables[8] = 1;
	variables[9] = 1;
}


void example2ODEs( concentrations &x , concentrations &dxdt , double t ){                                                         
    dxdt[0] = ((parameter[33] * x[7] * (x[9] - x[0]) / parameter[34]) / (1 + (x[9] / parameter[34]) + (x[0] / parameter[35])));
    dxdt[1] = -((parameter[27] * x[5] * (x[1] - x[8]) / parameter[28]) / (1 + (x[1] / parameter[28]) + (x[8] / parameter[29])));
    dxdt[2] = (parameter[0] / (1 + pow(x[0] / parameter[1],parameter[2]) + pow(parameter[3] / x[1],parameter[4]))) - (parameter[5] * x[2]);
    dxdt[3] = (parameter[6] / (1 + pow(x[0] / parameter[7],parameter[8]) + pow(parameter[9] / x[8],parameter[10]))) - (parameter[11] * x[3]);
    dxdt[4] = (parameter[12] / (1 + pow(x[0] / parameter[13],parameter[14]) + pow(parameter[15] / x[9],parameter[16]))) - (parameter[17] * x[4]);
    dxdt[5] = ((parameter[18] * x[2]) / (parameter[19] + x[2])) - (parameter[20] * x[5]);
    dxdt[6] = ((parameter[21] * x[3]) / (parameter[22] + x[3])) - (parameter[23] * x[6]);
    dxdt[7] = ((parameter[24] * x[3]) / (parameter[25] + x[4])) - (parameter[26] * x[7]);
	dxdt[8] = ((parameter[27] * x[5] * (x[1] - x[8]) / parameter[28]) / (1 + (x[1] / parameter[28]) + (x[8] / parameter[29]))) - ((parameter[30] * x[6] * (x[8] - x[9]) / parameter[31]) / (1 + (x[8] / parameter[31]) + (x[9] / parameter[32])));
	dxdt[9] = ((parameter[30] * x[6] * (x[8] - x[9]) / parameter[31]) / (1 + (x[8] / parameter[31]) + (x[9] / parameter[32]))) - ((parameter[33] * x[7] * (x[9] - x[0]) / parameter[34]) / (1 + (x[9] / parameter[34]) + (x[0] / parameter[35])));
}

void generateDataSet(double* parameters, concentrations &x, int varCount, double dt, double tf){
    double t = 0.0;
    int i;
    cout << t;
    for(i = 0; i < varCount; i++){
       cout << tab << x[i];
	}
	cout << endl;
    stepper_euler< concentrations > stepper;
    stepper.adjust_size( x );
    while(t <= tf){
        stepper.do_step( example2ODEs , x , t , dt );
		t += dt;
        cout << t;
        for(i = 0; i < varCount; i++){
	        cout << tab << x[i];
	     }
	     cout << endl;
    }
}

int main( int argc , char **argv ){
    const double dt = 0.1;
    concentrations x(VARIABLE_COUNT);
	initVariables(x);
	initParameters(parameter);
	generateDataSet(parameter,x,VARIABLE_COUNT,dt,10);
    return 0;
}
