#include <iostream>
#include <iterator>
#include <boost/numeric/odeint.hpp>
#include <boost/lambda/lambda.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

typedef vector< double > concentrations;

#define PARAMETER_COUNT 11
#define VARIABLE_COUNT 6
double parameter [PARAMETER_COUNT];

void initParameters(double* param, int column){
	switch(column){
		case 2:
	//With Constraint
	param[0] = 271;
	param[1] = 126;
	param[2] = 23;
	param[3] = 0.02;
	param[4] = 0.00005;
	param[5] = 0.00027;
	param[6] = 3.37;
	param[7] = 0.13;
	param[8] = 0.12;
	param[9] = 0.39;
	param[10] = 3;
	break;
		case 1:
	//Initial Velocity
	param[0] = 255;
	param[1] = 1511;
	param[2] = 145;
	param[3] = 2;
	param[4] = 0.5;
	param[5] = 3;
	param[6] = 0.16;
	param[7] = 0.46;
	param[8] = 0.04;
	param[9] = 0.3;
	param[10] = 1.9;
	break;
		default:
	//Without constraint
	param[0] = 204;
	param[1] = 142;
	param[2] = 18;
	param[3] = 0.35;
	param[4] = -0.5;
	param[5] = 0.0066;
	param[6] = -0.13;
	param[7] = 0.1;
	param[8] = 0.04;
	param[9] = 0.39;
	param[10] = 3.6;
	break;
	}
}

void initVariables(concentrations &variables){
	variables[0] = 5;
	variables[1] = 0;
	variables[2] = 0;
	variables[3] = 0;
	variables[4] = 1;
	variables[5] = 0;
}

double vGLK(concentrations &x, double* param){
	return (param[0] * x[0] * x[4]) / ((param[3] + x[0]) * (param[4] + x[4]));
}

double vPGI(concentrations &x, double* param){
	return (param[1] * (x[1] - (x[2] / param[9]))) / (x[1] + (param[5] * (1 + (x[2] / param[6]))));
}

double vPFKA(concentrations &x, double* param){
	return (param[2] * pow(x[2], param[10]) * x[4]) / ((pow(param[7], param[10]) + pow(x[2], param[10])) * (param[8] + x[4]));
}

void example1ODEs( concentrations &x , concentrations &dxdt , double t ){                                                         
    dxdt[0] = (-vGLK(x,parameter)) /100;
    dxdt[1] = (vGLK(x,parameter) - vPGI(x,parameter)) /100;
    dxdt[2] = (vPGI(x,parameter) - vPFKA(x,parameter)) /100;
    dxdt[3] = (vPFKA(x,parameter)) /100;
    dxdt[4] = (-vGLK(x,parameter) - vPFKA(x,parameter)) /100;
    dxdt[5] = (vGLK(x,parameter) + vPFKA(x,parameter)) /100;
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
        stepper.do_step( example1ODEs , x , t , dt );
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
	initParameters(parameter, 1);
	generateDataSet(parameter,x,VARIABLE_COUNT,dt,10);
    return 0;
}
