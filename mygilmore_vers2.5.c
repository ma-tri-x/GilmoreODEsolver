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
long double mu, BTait, sigma, omega, pv, pac, patm, C, c0, rho0v, phase, Rn1,Rn2,Rn;
long double rho0;

//------------------------------------------------------------------------------------------------
// asking whether you want to start with a volume 
// lower than that of the molecules
bool isPos (long double r){
    long double r3 = r*r*r;
    if (r3<=0.) return false;
    else return true;
}

// outsourced power law:
long double ThePower(long double r){
    long double r2 = r*r;
    long double r3 = r2*r;
    long double arg = (Rn*Rn*Rn)/(r3);
    if (arg<=0.) {
        cerr << "ERROR (ThePower): exiting because of nan-event. arg = " << arg << endl;
        cerr << "Rn = " << Rn << "  r = " << r << endl;
        exit(1);
    }
    return pow(arg,kappa);
}

// outsourced bubble pressure calculations:
vector<long double> pBubble (long double t, long double r, long double v){
    vector<long double> ReturnValues(3);
    long double pInfty = patm + pac*sin(omega*t + phase/180.*M_PI);
    long double dpInftyDt = pac*omega*cos(omega*t + phase/180.*M_PI);
    long double p = (patm + 2*sigma/Rn - pv)*ThePower(r) + pv -2*sigma/r - 4*mu*v/r;
    ReturnValues[0]=pInfty;
    ReturnValues[1]=dpInftyDt;
    ReturnValues[2]=p;
    return ReturnValues;
}

// the Gilmore equations: returning the acceleration
long double DGL (long double t, long double r, long double v)
{
    vector<long double> pValues(3);
    pValues = pBubble(t,r,v);
    long double pInfty = pValues[0];
    long double dpInftyDt = pValues[1];
    long double p = pValues[2];

    long double TaitArg = (p+BTait)/(pInfty+BTait);

    long double H = nTait/(nTait - 1.) * (pInfty+BTait)/rho0 * (pow(TaitArg,1.-1./nTait)-1.);

    // \dot H = H1*\dot p + H2*\dot pInfty
    long double H1 = pow(TaitArg,-1./nTait)/rho0;

    long double H2 =  1./(nTait - 1.) * (pow(TaitArg,1.-1./nTait) - nTait)/rho0;
    
    // \dot p = dp/dr * v + dp/d(dot r) * \ddot r; 
    // where dp/d(dot r) * \ddot r = 4mu/r \ddot r
    long double corrected_term = -3*kappa/r*v;// leftover of (derivation of ThePower wrt r) * v // see M.Sc. thesis p.7 but mind typo: R^3 is R^2.

    long double _dpdr_v = (patm + 2*sigma/Rn - pv)*ThePower(r)* corrected_term + (2.*sigma*v + 4.*mu*v*v)/r/r;
    //
    //
    C=sqrt(c0*c0 + (nTait-1.)*H);
    long double mach = v/C;

    long double FactorLHSto_ddotr = (1-mach)*(r+4.*mu/C*H1);

    long double a = ( -0.5*(3-mach)*v*v + (1+mach)*H + (1-mach)*r/C*  (    H1*_dpdr_v + H2*dpInftyDt    )  )    /   FactorLHSto_ddotr ;
    
    return a;
}

// RK method 3rd order with adaptive time step
void Runge_Kutta_embedded(long double &r, long double &v, long double &p, long double &dt, long double &t, long double epsilon) // functions introduced in "computational physics", Prof. Schmidt, Schleicher, SoSe2012
{
    long double old_r = r;
    long double old_v = v;
    
    long double k1v = DGL(t,  old_r,  old_v);
    long double k1r = old_v;
    
    long double k2v = DGL(t+dt,  old_r + dt*k1r,  old_v + dt*k1v);
    long double k2r = old_v + dt*k1v;
    
    long double k3v = DGL(t+dt/2.,  old_r + dt/4.*(k1r + k2r),  old_v + dt/4.*(k1v + k2v));
    long double k3r = old_v + dt/4.*(k1v+k2v);

    long double er = dt/3.* abs(k1r - 2.*k3r + k2r);
    long double ev = dt/3.* abs(k1v - 2.*k3v + k2v);
    
    // checking the errors and compare them to the relative
    // accuracy specified in main
    if (er < epsilon && ev < epsilon)
    {
        r = old_r + dt/6.*(k1r + 4.*k3r + k2r);
        v = old_v + dt/6.*(k1v + 4.*k3v + k2v);
        vector<long double> pValues(3);
        pValues = pBubble(t,r,v);
        p = pValues[2];
        t+=dt;
    }
    
    // time step is always updated
    // even if the "new" solution was
    // rejected
    if ( ev >= er )
    {
        dt = 0.9*dt*pow(epsilon/ev,1./3.);
    }else{
        dt = 0.9*dt*pow(epsilon/er,1./3.);
    }
}

int main(int argc, char *argv[])  // argumentscounter, charfeld, in das Eingaben gespeichert werden
{
    
    if (argc != 25 ) {               // im ersten Eintrag des char-Feldes steht der Programmname(?)
        cerr << "Usage: "<< argv[0] << " <Rstart> <vStart> <Rn> <deltaT> <Tstart> <Tend>" << endl;
        cerr                        << " <p_atm> <p_ac> <f> <mu> <sigma> " << endl;
        cerr                        << " <BTait> <nTait> <kappa> <pv> <rhoLiq> <SpecGasConst Bubble>" << endl;
        cerr                        << " <TempRef °C> <T0> <epsilon> <Rn2> <t1> <t2>" << endl;
        exit(1);
    }
    
    long double Rstart      = atof(argv[1]); // m
    long double RdotStart   = atof(argv[2]); // m/s
                Rn          = atof(argv[3]); // m
    long double deltaT      = atof(argv[4]); // s
    long double Tstart      = atof(argv[5]); // s
    long double Tend        = atof(argv[6]); // s
                patm        = atof(argv[7]); // Pa
                pac         = atof(argv[8]); // Pa
    long double f           = atof(argv[9]); // Hz
                    omega=2.*M_PI*f;
                mu          = atof(argv[10]); // Pa s
                sigma       = atof(argv[11]); // Pa m
                BTait       = atof(argv[12]); // Pa
                nTait       = atof(argv[13]); // dim.less
                kappa       = atof(argv[14]); // adiabat exp.
                pv          = atof(argv[15]); // Pa
                rho0        = atof(argv[16]); // kg/m^3
    long double SpecGasConst= atof(argv[17]); // J/(mol K) of vapour or gas at starting radius
    long double TempRef     = atof(argv[18]) + 273.15; // input in °C
    long double T0          = atof(argv[19]); // s;  scaling the time!
    long double epsilon     = atof(argv[20]);// was: 1e-2;
                phase       = atof(argv[21]);// was: 1e-2;
                Rn2         = atof(argv[22]); // m
    long double t1          = atof(argv[23]); // s
    long double t2          = atof(argv[24]); // s

    if (epsilon >=0.999999 || epsilon  <=0.){
        cerr << "epsilon is too high or negative";
        exit(1);
    }
    
    Rn1 = Rn;
    //phase = phase/180.*M_PI;   
    long double rho0vRef = patm/(SpecGasConst*TempRef); // kg/m^3 of vapour or gas at starting radius
                rho0v    = rho0vRef; // * (Rn1/Rstart)* (Rn1/Rstart)* (Rn1/Rstart);
    double GenGasConst = 8.3144621;
    cout << "rho_bubble,start = " << rho0v << " kg/m^3" << endl;
    cout << "speed of sound in bubble = " << sqrt(kappa*patm*pow(rho0v,kappa-1.)/pow(rho0vRef,kappa)) << "\n" << endl;

  
    /********************init!**********/
    long double r = Rstart;
    long double v = RdotStart;
    long double dt = deltaT;
    long double t = Tstart;
    c0 = sqrt((patm+BTait)*nTait/rho0); 
    long double p = pBubble(t,r,v)[2];
    cout << " initialized r,v,dt,c0 and p" << endl;
    /**************************************/
    cout << "time starts" << endl;
 
    ofstream plot;
    plot.open("gilmore2.dat");
    
    plot << "#   Rstart="       << Rstart       << "m" << endl;
    plot << "#   RdotStart="    << RdotStart    << "m/s" << endl;
    plot << "#   Rn="           << Rn           << "m" << endl;
    plot << "#   deltaT="       << deltaT       << "s" << endl;
    plot << "#   Tstart="       << Tstart       << "s" << endl;
    plot << "#   Tend="         << Tend         << "s" << endl;
    plot << "#   patm="         << patm         << "Pa" << endl;
    plot << "#   pac="          << pac          << "Pa" << endl;
    plot << "#   f="            << f            << "Hz" << endl;
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
 
    plot << "#t[s]\tr[m]\tv[m/s]\tdT [s]\tpBubble[Pa]\tp_ac[Pa]\tRn" << endl;
    plot << t << "\t" << r << "\t" << v << "\t" << dt << "\t" << ((patm + 2*sigma/Rn - pv)*ThePower(r)+pv) << "\t" << pBubble(t,r,v)[0] << "\t" << Rn << endl;

    while (t<Tend)
    {
        if (t >= t1 && t <= t2){
            Rn = (Rn2 - Rn1)/(t2 - t1)*(t - t1)+Rn1;
        }
        long double old_t = t;
        Runge_Kutta_embedded(r,v,p,dt,t,epsilon); // dt wird in der Funktion verändert
        if (isnan(r)||isnan(v)||isnan(p)||isnan(dt)|| isnan(t)){
            cerr << "ERROR-bla: exiting because of nan-event: r,v,p,dt,t:" << isnan(r)
                << isnan(v) << isnan(p) << isnan(dt) << isnan(t) << endl;
            break;
        }
        if(t!=old_t){
            plot << t << "\t" << r << "\t" << v << "\t" << dt << "\t" << ((patm + 2*sigma/Rn - pv)*ThePower(r)+pv) << "\t" << pBubble(t,r,v)[0] << "\t" << Rn << endl;
        }
    }
    plot.close();
}
