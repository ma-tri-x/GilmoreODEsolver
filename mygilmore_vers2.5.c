// #include "matrices.cpp"
#include <math.h>
// // #include "solvers.h"
// #include <iostream>
// #include <cstdlib>
// #include <string>
// #include <sstream>  // for things like ifstream but not with a file but with a string for input.
// #include <fstream>
// #include <stdio.h>  // für fstream
// #define GNUPLOT "gnuplot -persist"
using namespace std;

//############################################################################################//
//############# Program based on my 2012SoSe/Programmieren/sheet02 ###########################//
//############################################################################################//
//############# Runke Kutta 3rd order with adaptive time step handling #######################//
//############# solving the Gilmore Model with Van-der-Waals and acoustic pressure ###########//
//############################################################################################//
//############################################################################################//
//############################################################################################//
//############################################################################################//

#include <iostream>
#include <fstream> 
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
using namespace std;

// ----- global variables: -----------------------------
long double nTait, kappa;
long double sc_mu, sc_BTait, sc_sigma, sc_omega, U, sc_pv, sc_pac, sc_patm, sc_C, sc_c0, beta, sc_rho0v, phase;
const long double sc_rho0 = 1.;

//------------------------------------------------------------------------------------------------
// asking whether you want to start with a volume 
// lower than that of the molecules
bool isPos (long double r){
  long double r3 = r*r*r;
  if (r3-beta<=0.) return false;
  else return true;
}

// outsourced power law:
long double ThePower(long double r){
  long double r2 = r*r;
  long double r3 = r*r*r;
  long double arg = (1.-beta)/(r3-beta);
  if (arg<=0.) {
    cerr << "ERROR (ThePower): exiting because of nan-event.";
    exit(1);
  }
  return pow(arg,kappa);
}

// outsourced bubble pressure calculations:
vector<long double> pBubble (long double sc_t, long double r, long double v){
  vector<long double> ReturnValues(3);
  long double pInfty = sc_patm + sc_pac*sin(sc_omega*sc_t + phase/180.*M_PI);
  long double dpInftyDt = sc_pac*sc_omega*cos(sc_omega*sc_t + phase/180.*M_PI);
  long double sc_p = (sc_patm + 2*sc_sigma - sc_pv)*ThePower(r) + sc_pv -2*sc_sigma/r - 4*sc_mu*v/r;
//   cout << "das ganze ding=" << (sc_patm + 2*sc_sigma - sc_pv)*ThePower(r) << endl;
  ReturnValues[0]=pInfty;
  ReturnValues[1]=dpInftyDt;
  ReturnValues[2]=sc_p;
  return ReturnValues;
}

// the Gilmore equations: returning the acceleration
long double DGL (long double t, long double r, long double v)
{
  vector<long double> pValues(3);
  pValues = pBubble(t,r,v);
  long double pInfty = pValues[0];
  long double dpInftyDt = pValues[1];
  long double sc_p = pValues[2];

  long double TaitArg = (sc_p+sc_BTait)/(pInfty+sc_BTait);

  long double H = nTait/(nTait - 1.) * (pInfty+sc_BTait)/sc_rho0 * (pow(TaitArg,1.-1./nTait)-1.);

  long double H1 = pow(TaitArg,-1./nTait);//omitted: /sc_rho0

  long double H2 =  1./(nTait - 1.) * (pow(TaitArg,1.-1./nTait) - nTait);//omitted: /sc_rho0

  long double r2 = r*r;
  long double r3 = r*r*r;
  long double corrected_term = -3*kappa*r2*v/(r3-beta);// leftover of (derivation of ThePower wrt r) * v

  long double _dpdr_v = (sc_patm + 2*sc_sigma - sc_pv)*ThePower(r)* corrected_term + (2.*sc_sigma*v + 4.*sc_mu*v*v)/r/r;
  //
  //
  sc_C=sqrt(sc_c0*sc_c0 + (nTait-1.)*H);
  long double mach = v/sc_C;

  long double PseudoLHS = (1-mach)*(r+4.*sc_mu/sc_C*H1);

  long double a = -0.5*(3-mach)*v*v + (1+mach)*H + (1-mach)*r/sc_C*  (    H1*_dpdr_v + H2*dpInftyDt    );
  
  return a/PseudoLHS;
}

// RK method 3rd order with adaptive time step
void Runge_Kutta_embedded(long double &r, long double &v, long double &sc_p, long double &sc_dt, long double &sc_t, long double epsilon) // functions introduced in "computational physics", Prof. Schmidt, Schleicher, SoSe2012
{
 long double old_r = r;
 long double old_v = v;
 
 long double k1v = DGL(sc_t,  old_r,  old_v);
 long double k1r = old_v;
 
 long double k2v = DGL(sc_t+sc_dt,  old_r + sc_dt*k1r,  old_v + sc_dt*k1v);
 long double k2r = old_v + sc_dt*k1v;
 
 long double k3v = DGL(sc_t+sc_dt/2.,  old_r + sc_dt/4.*(k1r + k2r),  old_v + sc_dt/4.*(k1v + k2v));
 long double k3r = old_v + sc_dt/4.*(k1v+k2v);

 long double er = sc_dt/3.* abs(k1r - 2.*k3r + k2r);
 long double ev = sc_dt/3.* abs(k1v - 2.*k3v + k2v);
 
   // checking the errors and compare them to the relative
   // accuracy specified in main
 if (er < epsilon && ev < epsilon)
 {
   r = old_r + sc_dt/6.*(k1r + 4.*k3r + k2r);
   v = old_v + sc_dt/6.*(k1v + 4.*k3v + k2v);
   vector<long double> pValues(3);
   pValues = pBubble(sc_t,r,v);
   sc_p = pValues[2];
   sc_t+=sc_dt;
 }
 
   // time step is always updated
   // even if the "new" solution was
   // rejected
 if ( ev >= er )
 {
   sc_dt = 0.9*sc_dt*pow(epsilon/ev,1./3.);
 }else{
   sc_dt = 0.9*sc_dt*pow(epsilon/er,1./3.);
 }
}

// void Runge_Kutta_4th_order(long double &r, long double &v, long double &sc_p, long double sc_dt, long double &sc_t, long double epsilon) // etwas verwirrend hier: kjv sind ohne multiplikation mit dt (gemäß Vorlesung)
// {
//  long double old_r = r;
//  long double old_v = v;
//  
//  long double k1v = DGL(sc_t,old_r,old_v);
//  long double k1r = old_v;
//  
//  long double k2v = DGL(sc_t+sc_dt/2.,old_r+sc_dt/2.*k1r,old_v + k1v*sc_dt/2.);
//  long double k2r = old_v + k1v*sc_dt/2.;
//  
//  long double k3v = DGL(sc_t+sc_dt/2.,old_r + sc_dt/2.*k2r,old_v + sc_dt/2.*k2v);
//  long double k3r = old_v + sc_dt/2.*k2v;
//  
//  long double k4v = DGL(sc_t+sc_dt,old_r + sc_dt*k3r,old_v + sc_dt*k3v);
//  long double k4r = old_v + sc_dt*k3v;
//  
//    r = old_r + sc_dt/6.*(k1r + 2.*k2r + 2.*k3r + k4r);
//    v = old_v + sc_dt/6.*(k1v + 2.*k2v + 2.*k3v + k4v);
//    vector<long double> pValues(3);
//    pValues = pBubble(sc_t,r,v);
//    sc_p = pValues[2];
//    sc_t+=sc_dt;
//  }

int main(int argc, char *argv[])  // argumentscounter, charfeld, in das Eingaben gespeichert werden
{
  
   if (argc != 23 ) {               // im ersten Eintrag des char-Feldes steht der Programmname(?)
     cerr << "Usage: "<< argv[0] << " <Rstart> <vStart> <Rn> <deltaT> <Tstart> <Tend>" << endl;
     cerr                        << " <p_atm> <p_ac> <f> <bVan[m^3/kg]> <mu> <sigma> " << endl;
     cerr                        << " <BTait> <nTait> <kappa> <pv> <rhoLiq> <SpecGasConst Bubble>" << endl;
     cerr                        << " <TempRef °C> <T0> <epsilon> <Rn2> <t1> <t2>" << endl;
     exit(1);
   }else if (atof(argv[21])>=0.999999 || atof(argv[21])<=0.){
     cerr << "epsilon is too high or negative";
     exit(1);
   }
   
   long double Rstart      = atof(argv[1]); // m
   long double RdotStart   = atof(argv[2]); // m/s
   long double Rn          = atof(argv[3]); // m
   long double deltaT      = atof(argv[4]); // s
   long double Tstart      = atof(argv[5]); // s
   long double Tend        = atof(argv[6]); // s
   long double patm        = atof(argv[7]); // Pa
   long double pac         = atof(argv[8]); // Pa
   long double f           = atof(argv[9]); // Hz
      long double omega=2.*M_PI*f;
   long double bVan        = atof(argv[10]); // m^3/kg !!!!!
   long double mu          = atof(argv[11]); // Pa s
   long double sigma       = atof(argv[12]); // Pa m
   long double BTait       = atof(argv[13]); // Pa
   nTait                   = atof(argv[14]); // dim.less
   kappa                   = atof(argv[15]); // adiabat exp.
            cout << "kappa=" << kappa << endl;
   long double pv          = atof(argv[16]); // Pa
   long double rho0        = atof(argv[17]); // kg/m^3
   long double SpecGasConst= atof(argv[18]); // J/(mol K) of vapour or gas at starting radius
   long double TempRef     = 273.15 + atof(argv[19]); // input in °C
   long double T0          = atof(argv[20]); // s;  scaling the time!
   long double epsilon     = atof(argv[21]);// was: 1e-2;
               phase       = atof(argv[22]);// was: 1e-2;
   long double Rn2         = atof(argv[23]); // m
   long double t1          = atof(argv[24]); // s
   long double t2          = atof(argv[25]); // s
   cout << "phase=" << phase << endl;

   long double Rn1 = Rn;
       //phase = phase/180.*M_PI;   
       long double rho0vRef = patm/(SpecGasConst*TempRef); // kg/m^3 of vapour or gas at starting radius
       long double rho0v    = rho0vRef * (Rn/Rstart)* (Rn/Rstart)* (Rn/Rstart);
       double GenGasConst = 8.3144621;
       cout << "rho_bubble,start = " << rho0v << endl;
       cout << "speed of sound in bubble = " << sqrt(kappa*patm*pow(rho0v,kappa-1.)/pow(rho0vRef,kappa)) << "\n" << endl;

  
  /********************scaling!**********/
  cout << "scaling input data" << endl;
  long double r = Rstart/Rn;
  sc_omega = T0*omega;
  U = Rn/T0;
  long double v = RdotStart/U;
  long double sc_Tstart = Tstart/T0;
  long double sc_Tend = Tend/T0;
  long double sc_dt = deltaT/T0;
  long double pscale = rho0*U*U;
  sc_patm = patm/pscale;
  sc_pac  = pac/pscale;
  sc_pv   = pv/pscale;
  sc_BTait = BTait/pscale;
  sc_mu = mu/(pscale*Rn)*U;
  sc_sigma = sigma/(pscale*Rn);
  sc_c0 = sqrt((sc_patm+sc_BTait)*nTait); // bec. sc_rho0=1
  //   sc_rho0v = rho0v/rho0; // not needed
  beta = bVan*patm/(GenGasConst*TempRef); // beta stays beta because unitless
  //beta = 0.0016;
  /**************************************/
  
  /**************** asking for negativity Rstart ********/
    if (!isPos(r)){
      cerr << "ERROR: Startvolume < Hard Cores' Volume" << endl;
      exit(1);
    }
  /******************************************************/
 long double sc_t = sc_Tstart;
 long double sc_p = pBubble(sc_t,r,v)[2];
 
 /*long double epsilon=1e-2;*/ // accuracy of time step adapting
 
 cout << "time starts" << endl;
 cout << "(V_start-V_Hardcores)/(4/3*pi) = " << Rstart*Rstart*Rstart*(1.-beta) << endl;
 
 ofstream plot;
//  ofstream plot2;
 plot.open("gilmore2.dat");
//  plot2.open("gilmore_unitless.dat");
 plot << "#   Rstart="       << Rstart       << "m" << endl;
 plot << "#   RdotStart="    << RdotStart    << "m/s" << endl;
 plot << "#   Rn="           << Rn           << "m" << endl;
 plot << "#   deltaT="       << deltaT       << "s" << endl;
 plot << "#   Tstart="       << Tstart       << "s" << endl;
 plot << "#   Tend="         << Tend         << "s" << endl;
 plot << "#   patm="         << patm         << "Pa" << endl;
 plot << "#   pac="          << pac          << "Pa" << endl;
 plot << "#   f="            << f            << "Hz" << endl;
 plot << "#   bVan="         << bVan         << "m³/kg" << endl;
 plot << "#   mu="           << mu           << "Pa s" << endl;
 plot << "#   sigma="        << sigma        << "Pa m (I think)" << endl;
 plot << "#   BTait="        << BTait        << "Pa" << endl;
 plot << "#   nTait="        << nTait        << "no unit" << endl;
 plot << "#   kappa="        << kappa        <<  endl;
 plot << "#   pv="           << pv           << "Pa" << endl;
 plot << "#   rho0="         << rho0         << "kg/m³" << endl;
 plot << "#   SpecGasConst=" << SpecGasConst << "J/molK" << endl;
 plot << "#   TempRef="      << TempRef      << "deg Celsius." << endl;
 plot << "#   T0="           << T0           <<  endl;
 plot << "#   Rn2"           << Rn2          << "m"; // m
 plot << "#   t1"            << t1           << "s"; // s
 plot << "#   t2"            << t2           << "s"; // s
 plot << "# "                                << endl;
 
 plot << "#t[s]\tr[m]\tv[m/s]\tdT [s]\tpBubble[Pa]\tp_ac[Pa]" << endl;
 plot << sc_t*T0 << "\t" << r*Rn << "\t" << v*U << "\t" << sc_dt*T0 << "\t" << ((sc_patm + 2*sc_sigma - sc_pv)*ThePower(r)+sc_pv)*pscale << "\t" << pBubble(sc_t,r,v)[0]*pscale << endl;
//  plot2 << "#t\tr\tv\tdt\tpBubble" << endl;
//  plot2 << sc_t << "\t" << r << "\t" << v << "\t" << sc_dt << "\t" << sc_p << endl;
 while (sc_t<sc_Tend)
   {
    if (sc_t*T0 >= t1 && sc_t*T0 <= t2){
        Rn = (Rn2 - Rn1)/(t2 - t1)*(sc_t*T0 - t1);
        /********************re-scaling!**********/
        r = r*Rn1/Rn;
        U = Rn/T0;
        v = (v*Rn1/T0)/U;
        pscale = rho0*U*U;
        sc_patm = patm/pscale;
        sc_pac  = pac/pscale;
        sc_pv   = pv/pscale;
        sc_BTait = BTait/pscale;
        sc_mu = mu/(pscale*Rn)*U;
        sc_sigma = sigma/(pscale*Rn);
        sc_c0 = sqrt((sc_patm+sc_BTait)*nTait); // bec. sc_rho0=1
        /**************************************/
    }
	long double old_t = sc_t;
	  Runge_Kutta_embedded(r,v,sc_p,sc_dt,sc_t,epsilon); // dt wird in der Funktion verändert
	if (isnan(r)||isnan(v)||isnan(sc_p)||isnan(sc_dt)|| isnan(sc_t)){
	  cerr << "ERROR: exiting because of nan-event" << endl;
	  break;
	}
	if(sc_t!=old_t){
	  plot << sc_t*T0 << "\t" << r*Rn << "\t" << v*U << "\t" << sc_dt*T0 << "\t" << ((sc_patm + 2*sc_sigma - sc_pv)*ThePower(r)+sc_pv)*pscale << "\t" << pBubble(sc_t,r,v)[0]*pscale << endl;
// 	  plot2 << sc_t << "\t" << r << "\t" << v << "\t" << sc_dt << "\t" << sc_p << endl;
	}
   }
 plot.close();
//  plot2.close();
}
