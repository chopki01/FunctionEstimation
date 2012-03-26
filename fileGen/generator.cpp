#include <iostream>
#include <iterator>
#include <boost/numeric/odeint.hpp>
#include <boost/lambda/lambda.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

typedef vector< double > concentrations;

int met = 0;
const int GLUC = met++;
const int GLUC6P = met++;
const int ATP = met++;
const int ADP = met++;
const int PGL = met++;
const int NAD = met++;
const int NADH = met++;
const int NADP = met++;
const int NADPH = met++;
const int PG = met++;
const int KDPG = met++;
const int G3P = met++;
const int GAP = met++;
const int PYR = met++;
const int BPG = met++;
const int G2P = met++;
const int PEP = met++;
const int ACET = met++;
const int ETOH = met++;

const int VARIABLE_COUNT = met;

// 1 = Glucokinase;
const double kcat_1 = 600;
const double Km_GLUC_1 = 0.22;
const double Km_ATP_1 = 0.8;
const double Ki_GLUC6P_1 = 15;
// 2 = Glucose-6-P;
const double kcat_NADP_2 = 292.5;
const double kcat_NAD_2 = 510;
const double Km_GLUC6P_2 = 0.17;
const double Km_NADP_2 = 0.04;
const double Km_NAD_2 = 0.21;
const double Ki_ATP_2 = 1.4;
// 3 = 6-phosphogluconolactonase;
const double kcat_3 = 6600;
const double Km_PGL_3 = 0.025;
const double Ki_GLUC6P_3 = 0.3;
// 4 = 6-phosphogluconate;
const double kcat_4 = 367.5;
const double Km_PG_4 = 0.04;
const double Ki_G3P_4 = 2;
// 5 = 2-keto-3-deoxy-6-phosphogluconate;
const double kcat_5 = 900;
const double Km_KDPG_5 = 0.25;
// 6 = Glyceraldehyde-3-P;
const double kcat_f_6 = 25;
const double kcat_r_6 = 307.5;
const double Km_GAP_6 = 1;
const double Km_NAD_6 = 0.09;
const double Km_BPG_6 = 1;
const double Km_NADH_6 = 0.02;
// 7 = 3-phosphoglycerate;
const double kcat_f_7 = 8.4;
const double kcat_r_7 = 1200;
const double Km_ADP_7 = 0.8;
const double Km_BPG_7 = 1;
const double Km_ATP_7 = 0.8;
const double Km_G3P_7 = 1;
// 8 = Phosphoglycerate;
const double kcat_8 = 3000;
const double Km_G3P_8 = 1.1;
// 9 = Enolase;
const double kcat_9 = 240;
const double Km_G2P_9 = 0.08;
// 10 = Pyruvate;
const double kcat_10 = 450;
const double Km_PEP_10 = 0.08;
const double Km_ADP_10 = 0.17;
// 11 = Pyruvate;
const double kcat_11 = 181;
const double Km_PYR_11 = 0.4;
// 12 = Alcohol;
const double kcat_f_12 = 2300;
const double kcat_r_12 = 190;
const double Km_ACET_12 = 0.086;
const double Km_NADH_12 = 0.027;
const double Ki_NADH_12 = 0.0076;
const double Km_ETOH_12 = 4.8;
const double Km_NAD_12 = 0.073;
const double Ki_NAD_12 = 0.024;

void initVariables(concentrations &variables){
	for(int i = 0; i < VARIABLE_COUNT; i++){
		variables[i] = 0;
	}
	variables[GLUC] = 2;
	variables[ATP] = 2;
	variables[ADP] = 1;
	variables[NAD] = 0;
	variables[NADH] = 0;
	variables[NADP] = 0;
	variables[NADPH] = 0;
	variables[ACET] = 0;
	variables[ETOH] = 2;
}

double vGK(concentrations &x){
	double n = (x[GLUC] / Km_GLUC_1) * (x[ATP] / (Km_ATP_1 * (1 + (x[GLUC6P] / Ki_GLUC6P_1))));
	double d = 1 + (x[GLUC] / Km_GLUC_1);
	d += x[ATP] / (Km_ATP_1 * (1 + (x[GLUC6P] / Ki_GLUC6P_1)));
	d += (x[GLUC] * x[ATP]) / (Km_GLUC_1 * Km_ATP_1 * (1 + (x[GLUC6P] / Ki_GLUC6P_1)));
	d += (x[GLUC6P] / Ki_GLUC6P_1);
	return n /d;
}

double vGPD(concentrations &x){
	double gluc6p = x[GLUC6P] / Km_GLUC6P_2;
	double nadp = x[NADP] / Km_NADP_2;
	double atp = x[ATP] / Ki_ATP_2;
	double nad = x[NAD] / (Km_NAD_2 * (1 + atp));
	double result = (gluc6p * nadp) / (1 + gluc6p + nadp + gluc6p * nadp);
	result += (gluc6p * nad) / (1 + gluc6p + nad + gluc6p * nad + atp);
	return result;
}

double vPGLS(concentrations &x){
	double gluc6p = x[GLUC6P] / Ki_GLUC6P_3;
	double pgl = x[PGL] / (Km_PGL_3 * (1 + gluc6p));
	return pgl / (1 + pgl);
}

double vPGD(concentrations &x){
	double g3p = x[G3P] / Ki_G3P_4;
	double pg = x[PG] / (Km_PG_4 * (1 + g3p));
	return pg / (1 + pg);
}

double vKDPGA(concentrations &x){
	double kdpg = x[KDPG] / Km_KDPG_5;
	return kdpg / (1 + kdpg);
}

double vGAPD(concentrations &x){
	double gap = x[GAP] / Km_GAP_6;
	double nad = x[NAD] / Km_NAD_6;
	double bpg = x[BPG] / Km_BPG_6;
	double nadh = x[NADH] / Km_NADH_6;
	return (0.0813 * gap * nad - bpg * nadh) / (1 + gap + nad + gap * nad + bpg + nadh + bpg * nadh);
}

double vG3PK(concentrations &x){
	double adp = x[ADP] / Km_ADP_7;
	double bpg = x[BPG] / Km_BPG_7;
	double atp = x[ATP] / Km_ATP_7;
	double g3p = x[G3P] / Km_G3P_7;
	return (0.007 * adp * bpg - atp * g3p) / (1 + adp + bpg + adp * bpg + atp + g3p + atp * g3p);
}

double vGPM(concentrations &x){
	double g3p = x[G3P] / Km_G3P_8;
	return g3p / (1 + g3p);
}

double vENO(concentrations &x){
	double g2p = x[G2P] / Km_G2P_9;
	return g2p / (1 + g2p);
}

double vPYRK(concentrations &x){
	double pep = x[PEP] / Km_PEP_10;
	double adp = x[ADP] / Km_ADP_10;
	return (pep * adp) / (1 + pep + adp + pep * adp);
}

double vPYRD(concentrations &x){
	double pyr = x[PYR] / Km_PYR_11;
	return pyr / (1 + pyr);
}

double vADHI(concentrations &x){
	double acet = x[ACET] / Km_ACET_12;
	double nadh = x[NADH] / Km_NADH_12;
	double nadhi = x[NADH] / Ki_NADH_12;
	double etoh = x[ETOH] / Km_ETOH_12;
	double nad = x[NAD] / Km_NAD_12;
	double nadi = x[NAD] / Ki_NAD_12;
	return (acet * nadh - 0.0826 * etoh * nad) / (1 + acet + nadh + acet * nadhi + etoh + nad + etoh * nadi);
}

void example3ODEs(concentrations &x , concentrations &dxdt , double t ){
	double vGLUC = 0;
	double vETOH = 0;
	double v1 = vGK(x);
	double v2 = vGPD(x);
	double v3 = vPGLS(x);
	double v4 = vPGD(x);
	double v5 = vKDPGA(x);
	double v6 = vGAPD(x);
	double v7 = vG3PK(x);
	double v8 = vGPM(x);
	double v9 = vENO(x);
	double v10 = vPYRK(x);
	double v11 = vPYRD(x);
	double v12 = vADHI(x);
	for(int i = 0; i < VARIABLE_COUNT; i++){
		dxdt[i] = 0;
	}
	dxdt[GLUC] = vGLUC - v1;
	dxdt[GLUC6P] = v1 - v2;
	dxdt[PGL] = v2 - v3;
	dxdt[PG] = v3 - v4;
	dxdt[KDPG] = v4 - v5;
	dxdt[GAP] = v5 - v6;
	dxdt[BPG] = v6 - v7;
	dxdt[G3P] = v7 - v8;
	dxdt[G2P] = v8 - v9;
	dxdt[PEP] = v9 - v10;
	dxdt[PYR] = v5 + v10 - v11;
	dxdt[ACET] = v11 - v12;
	dxdt[ETOH] = v12 - vETOH;
}

void generateDataSet(concentrations &x, int varCount, double dt, double tf){
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
		stepper.do_step( example3ODEs , x , t , dt );
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
	generateDataSet(x,VARIABLE_COUNT,dt,10);
	return 0;
}
