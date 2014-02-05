/**
   Upscale parallell model analytically in the capillary/viscous limit
 */

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
//#include <ctime>
#include <cmath>
#include <cfloat> // for DBL_MAX/DBL_MIN
#include <map>
#include <sys/utsname.h>

#include <opm/core/utility/MonotCubicInterpolator.hpp>

using namespace Opm;
using namespace std;


double arithmeticAvg(vector<double> parameter, vector<double> PV) {
    double avg = (PV[0]*parameter[0] + PV[1]*parameter[1]) / (PV[0] + PV[1]);
    return avg;
}

double harmonicAvg(vector<double> parameter, vector<double> PV) {
    double avg = (PV[0] + PV[1]) / (PV[0]/parameter[0] + PV[1]/parameter[1]);
    return avg;
}

vector<double> calcSatCE(double sat_init, vector<double> PV, vector<double> Jfactor,
			 MonotCubicInterpolator Jfunctions, MonotCubicInterpolator InvJfunctions) {
    vector<double> sat(2, 0.0);
    double J = Jfunctions.evaluate(sat_init);
    vector<double> Pc;
    Pc.push_back(J/Jfactor[0]);
    Pc.push_back(J/Jfactor[1]);
    double upscaled_Pc = arithmeticAvg(Pc, PV);
    vector<double> Jvec;
    Jvec.push_back(Jfactor[0]*upscaled_Pc);
    Jvec.push_back(Jfactor[1]*upscaled_Pc);
    sat[0] = InvJfunctions.evaluate(Jvec[0]);
    sat[1] = InvJfunctions.evaluate(Jvec[1]);
    return sat;
}

vector<MonotCubicInterpolator> calcFracFlowFunctions(vector<MonotCubicInterpolator> KrfunctionsW, 
						     vector<MonotCubicInterpolator> KrfunctionsO, 
						     vector<double> visc, int points) {
    // Create v_w/v_o as function of S
    vector<double> swir, swor; 
    swir.push_back(KrfunctionsW[0].getMinimumX().first);
    swir.push_back(KrfunctionsW[1].getMinimumX().first);
    swor.push_back(KrfunctionsO[0].getMaximumX().first);
    swor.push_back(KrfunctionsO[1].getMaximumX().first);

    vector<double> sat1(points, 0.0), sat2(points, 0.0);
    vector<double> relpermratio1(points, 0.0);
    vector<double> relpermratio2(points, 0.0);
    for (int i = 0; i < points; ++i) {
	sat1[i] = swir[0] + i*(swor[0]-swir[0])/(points-1);
	sat2[i] = swir[1] + i*(swor[1]-swir[1])/(points-1);
	relpermratio1[i] = KrfunctionsW[0].evaluate(sat1[i]) / KrfunctionsO[0].evaluate(sat1[i]);
	relpermratio2[i] = KrfunctionsW[1].evaluate(sat2[i]) / KrfunctionsO[1].evaluate(sat2[i]);
    }
    
    vector<MonotCubicInterpolator> FracFlowFun;
    FracFlowFun.push_back(MonotCubicInterpolator(sat1, relpermratio1));
    FracFlowFun.push_back(MonotCubicInterpolator(sat2, relpermratio2));
    FracFlowFun[0].scaleData(visc[1]/visc[0]);
    FracFlowFun[1].scaleData(visc[1]/visc[0]);
    
    if (FracFlowFun[0].isStrictlyMonotone()) {
	FracFlowFun.push_back(MonotCubicInterpolator(FracFlowFun[0].get_fVector(), FracFlowFun[0].get_xVector()));
    }
    else {
	cerr << "Frac flow function for rock 1 not invertible!" << endl;
    }
    if (FracFlowFun[1].isStrictlyMonotone()) {
	FracFlowFun.push_back(MonotCubicInterpolator(FracFlowFun[1].get_fVector(), FracFlowFun[1].get_xVector()));
    }
    else {
	cerr << "Frac flow function for rock 2 not invertible!" << endl;
    }    
    
    return FracFlowFun;
}


vector<double> calcSatVL(double sat_init, vector<double> PV, vector<double> visc, 
			 vector<MonotCubicInterpolator> FracFlowFun, vector<MonotCubicInterpolator> InvFracFlowFun) {
    vector<double> sat(2, 0.0);
    vector<double> fracflow;

    vector<double> swir, swor; 
    swir.push_back(FracFlowFun[0].getMinimumX().first);
    swir.push_back(FracFlowFun[1].getMinimumX().first);
    swor.push_back(FracFlowFun[0].getMaximumX().first);
    swor.push_back(FracFlowFun[1].getMaximumX().first);  

    fracflow.push_back(FracFlowFun[0].evaluate(sat_init));
    fracflow.push_back(FracFlowFun[1].evaluate(sat_init));
    double fracflowpoint = arithmeticAvg(fracflow, PV);

    if (sat_init == swor[0])       sat[0] = swor[0];
    else if (sat_init == swir[0])  sat[0] = swir[0];
    else                           sat[0] = InvFracFlowFun[0].evaluate(fracflowpoint);
   
    if (sat_init == swor[1])       sat[1] = swor[1];
    else if (sat_init == swir[1])  sat[1] = swir[1];
    else                           sat[1] = InvFracFlowFun[1].evaluate(fracflowpoint);

    return sat;
}


int main(int varnum, char** vararg)
{

    /*
      Populate options-map with default values
    */
    map<string,string> options;
    options.insert(make_pair("points",             "31"   )); // Number of saturation points (uniformly distributed within saturation endpoints)
    options.insert(make_pair("relPermCurve",       "2")); // Which column in the rock types are upscaled
    options.insert(make_pair("jFunctionCurve",     "4")); // Which column in the rock type file is the J-function curve
    options.insert(make_pair("surfaceTension",     "11")); // Surface tension given in dynes/cm
    options.insert(make_pair("output",             "")); // If this is set, output goes to screen and to this file. 
    // options.insert(make_pair("gravity",            "0.0")); // default is no gravitational effects
    // options.insert(make_pair("waterDensity",       "1.0")); // default density of water, only applicable to gravity
    // options.insert(make_pair("oilDensity",         "0.6")); // ditto
    options.insert(make_pair("outputprecision",    "10")); // number of significant numbers to print
    options.insert(make_pair("maxPermContrast",    "1e7")); // maximum allowed contrast in each single-phase computation
    options.insert(make_pair("minPerm",            "1e-12")); // absolute minimum for allowed cell permeability
    options.insert(make_pair("maxPerm",            "100000")); // maximal allowed cell permeability
    options.insert(make_pair("minPoro",            "0.0001")); // this limit is necessary for pcmin/max computation
    options.insert(make_pair("saturationThreshold","0.00001")); // accuracy threshold for saturation, we ignore Pc values that
    options.insert(make_pair("perm1",              "100")); // Permeability of rock 1
    options.insert(make_pair("perm2",              "1")); // Permeability of rock 2
    options.insert(make_pair("poro1",              "0.1")); // Porosity of rock 1
    options.insert(make_pair("poro2",              "0.1")); // Porosity of rock 2
    options.insert(make_pair("viscW",              "3e-4")); // Viscosity water
    options.insert(make_pair("viscO",              "2e-3")); // Viscosity oil
    options.insert(make_pair("useRockFiles",       "false")); // Use rock files. If not Brooks-Corey type Cp and Corey type relperm is used
    options.insert(make_pair("rock1File",          "")); // Rock 1 curves
    options.insert(make_pair("rock2File",          "")); // Rock 2 curves
    options.insert(make_pair("volFrac1",           "0.5")); // Volume fraction of rock 1

    // Conversion factor, multiply mD numbers with this to get m² numbers
    const double milliDarcyToSqMetre = 9.869233e-16;
    // Reference: http://www.spe.org/spe-site/spe/spe/papers/authors/Metric_Standard.pdf

    /* Loop over all command line options in order to look 
       for options. 

       argidx loops over all the arguments here, and updates the
       variable 'argeclindex' *if* it finds any legal options,
       'argeclindex' is so that vararg[argeclindex] = the eclipse
       filename. If options are illegal, argeclindex will be wrong, 
      
    */
    for (int argidx = 1; argidx < varnum; argidx += 2)  {
	if (string(vararg[argidx]).substr(0,1) == "-")    {
	    string searchfor = string(vararg[argidx]).substr(1); // Chop off leading '-'
	    /* Check if it is a match */
	    if (options.count(searchfor) == 1) {
		options[searchfor] = string(vararg[argidx+1]);
		cout << "Parsed command line option: " << searchfor << " := " << vararg[argidx+1] << endl;
	    }
	    else {
		cout << "Option -" << searchfor << " unrecognized." << endl;
		//usageandexit();
	    }
	}
    }

    /* Create/load rock parameters */
    std::vector<MonotCubicInterpolator> Jfunctions; // Holds the loaded J-functions.
    std::vector<MonotCubicInterpolator> InvJfunctions; // Holds the inverse of the loaded J-functions.
    std::vector<MonotCubicInterpolator> KrfunctionsW; // Holds relperm-curves for phase 1 for water
    std::vector<MonotCubicInterpolator> KrfunctionsO; // Holds relperm-curves for phase 2 for oil
    std::vector<string> JfunctionNames; // Placeholder for the names of the loaded J-functions.

    if (options["useRockFiles"] == "true") { // Use input files
	std::vector<string> ROCKFILENAMES;
	ROCKFILENAMES.push_back(options["rock1File"]);
	ROCKFILENAMES.push_back(options["rock2File"]);
	const int jFunctionCurve        = atoi(options["jFunctionCurve"].c_str());
	const int relPermCurve = atoi(options["relPermCurve"].c_str());

	for (int rock = 0; rock < 2; ++rock) { // Loop trough rock files
	    const string rockFile = ROCKFILENAMES[rock];

	    // Check if rock file is readable
	    ifstream rockstream(rockFile, ios::in);
	    if (rockstream.fail()) {
		cerr << "Error: Filename " << rockFile << " not found or not readable." << endl;
	    }
	    rockstream.close(); 
	
	    MonotCubicInterpolator Jtmp;
	    try {
		Jtmp = MonotCubicInterpolator(rockFile, 1, jFunctionCurve); 
	    }
	    catch (const char * errormessage) {
		cerr << "Error: " << errormessage << endl;
		cerr << "Check filename and -jFunctionCurve" << endl;
	    }     
	    // Invert J-function, now we get saturation as a function of pressure:
	    if (Jtmp.isStrictlyMonotone()) {
		Jfunctions.push_back(Jtmp);
		InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
		JfunctionNames.push_back(rockFile);
		KrfunctionsW.push_back(MonotCubicInterpolator(rockFile, 1, relPermCurve));
		KrfunctionsO.push_back(MonotCubicInterpolator(rockFile, 1, relPermCurve+1));
	    }
	    else {
		cerr << "Error: Jfunction " << rock+1 << " in rock file " << rockFile << " was not invertible." << endl;
	    }
	}
    }
    else { // Use Cp/relperm from closed form
	const int p = 101; // Points to include
	vector<double> s(p,0.0);
	vector<double> krw(p,0.0);
	vector<double> kro(p,0.0);
	vector<double> Jfun(p,0.0);
	for (int i = 0; i < p; ++i) {
	    s[i]    = i/(double(p-1.0));
	    krw[i]  = s[i]*s[i];
	    kro[i]  = (1-s[i])*(1-s[i]);
	    if (s[i] == 0.0) { Jfun[i] = sqrt(DBL_MAX); }
	    else if (s[i] == 1.0) { Jfun[i] = -sqrt(DBL_MAX); }
	    else {Jfun[i] = 1/(s[i]*s[i]) - 1/((1-s[i])*(1-s[i])); }
	}
	KrfunctionsW.push_back(MonotCubicInterpolator(s, krw));
	KrfunctionsW.push_back(MonotCubicInterpolator(s, krw));
	KrfunctionsO.push_back(MonotCubicInterpolator(s, kro));
	KrfunctionsO.push_back(MonotCubicInterpolator(s, kro)); // Use same curves for both
	MonotCubicInterpolator Jtmp(s, Jfun); 
	// Invert J-function, now we get saturation as a function of pressure:
	if (Jtmp.isStrictlyMonotone()) {
	    Jfunctions.push_back(Jtmp);
	    Jfunctions.push_back(Jtmp);
	    InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
	    InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
	    JfunctionNames.push_back("Brooks-Corey");
	    JfunctionNames.push_back("Brooks-Corey");
	}
	else {
	    cerr << "Error: Jfunction (closed form) was not invertible." << endl;
	}
    }
    
    /* Pre-calculations */
    double st = atof(options["surfaceTension"].c_str())*0.001;
    vector<double> perm, poro, visc, PV, Jfactor, volFrac;
    perm.push_back(atof(options["perm1"].c_str()));
    perm.push_back(atof(options["perm2"].c_str()));
    poro.push_back(atof(options["poro1"].c_str()));
    poro.push_back(atof(options["poro2"].c_str()));
    visc.push_back(atof(options["viscW"].c_str()));
    visc.push_back(atof(options["viscO"].c_str()));
    volFrac.push_back(atof(options["volFrac1"].c_str()));
    volFrac.push_back(1 - volFrac[0]);
    PV.push_back(poro[0]*volFrac[0]);
    PV.push_back(poro[1]*volFrac[1]);
    Jfactor.push_back(sqrt(milliDarcyToSqMetre*perm[0]/poro[0])/st);
    Jfactor.push_back(sqrt(milliDarcyToSqMetre*perm[1]/poro[1])/st);
    vector<MonotCubicInterpolator> FracFlowFunTmp = calcFracFlowFunctions(KrfunctionsW, KrfunctionsO, visc, 1000);
    vector<MonotCubicInterpolator> FracFlowFun, InvFracFlowFun;
    FracFlowFun.push_back(FracFlowFunTmp[0]);
    FracFlowFun.push_back(FracFlowFunTmp[1]);
    InvFracFlowFun.push_back(FracFlowFunTmp[2]);
    InvFracFlowFun.push_back(FracFlowFunTmp[3]);

    stringstream outputtmp;

    /* Upscale abs perm */
    double upscaled_perm_a = arithmeticAvg(perm, PV);
    double upscaled_perm_h = harmonicAvg(perm, PV);

    /* Upscale rel perm */
    const int points = atoi(options["points"].c_str());
    vector<double> zeros(points, 0.0);
    vector<double> upscaled_sat_CE = zeros;
    vector<double> upscaled_sat_VL = zeros;
    vector<double> upscaled_krw_CE_a = zeros;
    vector<double> upscaled_kro_CE_a = zeros;
    vector<double> upscaled_krw_CE_h = zeros;
    vector<double> upscaled_kro_CE_h = zeros;
    vector<double> upscaled_krw_VL_a = zeros;
    vector<double> upscaled_kro_VL_a = zeros;
    vector<double> upscaled_krw_VL_h = zeros;
    vector<double> upscaled_kro_VL_h = zeros;

    for (int i = 0; i < points; ++i) {
	double sat_init = i/(double(points-1.0));
	vector<double> kw(2, 0.0), ko(2, 0.0);
	
	// CE
	vector<double> sat_CE = calcSatCE(sat_init, PV, Jfactor, Jfunctions[0], InvJfunctions[0]);
	upscaled_sat_CE[i] = arithmeticAvg(sat_CE, PV);
	kw[0] = KrfunctionsW[0].evaluate(sat_CE[0])*perm[0];
	kw[1] = KrfunctionsW[1].evaluate(sat_CE[1])*perm[1];
	ko[0] = KrfunctionsO[0].evaluate(sat_CE[0])*perm[0];
	ko[1] = KrfunctionsO[1].evaluate(sat_CE[1])*perm[1];
	upscaled_krw_CE_a[i] = arithmeticAvg(kw, PV)/upscaled_perm_a;
	upscaled_kro_CE_a[i] = arithmeticAvg(ko, PV)/upscaled_perm_a;
	upscaled_krw_CE_h[i] = harmonicAvg(kw, PV)/upscaled_perm_h;
	upscaled_kro_CE_h[i] = harmonicAvg(ko, PV)/upscaled_perm_h;
	
	// VL
	vector<double> sat_VL = calcSatVL(sat_init, PV, visc, FracFlowFun, InvFracFlowFun);
	upscaled_sat_VL[i] = arithmeticAvg(sat_VL, PV);
	kw[0] = KrfunctionsW[0].evaluate(sat_VL[0])*perm[0];
	kw[1] = KrfunctionsW[1].evaluate(sat_VL[1])*perm[1];
	ko[0] = KrfunctionsO[0].evaluate(sat_VL[0])*perm[0];
	ko[1] = KrfunctionsO[1].evaluate(sat_VL[1])*perm[1];
	upscaled_krw_VL_a[i] = arithmeticAvg(kw, PV)/upscaled_perm_a;
	upscaled_kro_VL_a[i] = arithmeticAvg(ko, PV)/upscaled_perm_a;
	upscaled_krw_VL_h[i] = harmonicAvg(kw, PV)/upscaled_perm_h;
	upscaled_kro_VL_h[i] = harmonicAvg(ko, PV)/upscaled_perm_h;
    }
  
    outputtmp << "######################################################################" << endl
	      << "# Analytical steadystate upscaling (Arithmetic and Harmonic average)" << endl
	      << "######################################################################" << endl
	      << "# Parallel model with two alternating rocks:" << endl
	      << "#   Rock 1: perm = " << perm[0] << "mD, poro = " << poro[0] << ", volume fraction = " << volFrac[0]
	      << ", curves = \'" << JfunctionNames[0] << "\' (" << Jfunctions[0].getSize() << " points)" << endl 
	      << "#   Rock 2: perm = " << perm[1] << "mD, poro = " << poro[1] << ", volume fraction = " << volFrac[1]
	      << ", curves = \'" << JfunctionNames[1] << "\' (" << Jfunctions[1].getSize() << " points)" << endl
	      << "######################################################################" << endl
	      << "# Upscaled absolute permeability (arithmetic): " << upscaled_perm_a << "mD" << endl
	      << "# Upscaled absolute permeability (harmonic):   " << upscaled_perm_h << "mD" << endl
	      << "######################################################################" << endl
	      << "# Capillary limit\t\t\t\t\t\t\t\tViscous limit" << endl
	      << "# Sw\t\tkrw_a\t\tkro_a\t\tkrw_h\t\tkro_h\t\tSw\t\tkrw_a\t\tkro_a\t\tkrw_h\t\tkro_h" << endl;
    outputtmp << setprecision(atoi(options["outputprecision"].c_str())) << fixed;
    for (int i = 0; i < points; ++i) {
	outputtmp << upscaled_sat_CE[i] << "\t" << upscaled_krw_CE_a[i] << "\t" << upscaled_kro_CE_a[i] << "\t" 
		  << upscaled_krw_CE_h[i] << "\t" << upscaled_kro_CE_h[i] << "\t" 
		  << upscaled_sat_VL[i] << "\t" << upscaled_krw_VL_a[i] << "\t" << upscaled_kro_VL_a[i] << "\t" 
		  << upscaled_krw_VL_h[i] << "\t" << upscaled_kro_VL_h[i] << endl;
    }

    /* Print to screen (and file)*/
    cout << outputtmp.str();

    if (options["output"] != "") {
	cout << "######################################################################" << endl
	     << "Writing results to " << options["output"] << endl;
	ofstream outfile;
	outfile.open(options["output"].c_str(), ios::out | ios::trunc);
	outfile << outputtmp.str();
	outfile.close();
    }

    return 0;
}

