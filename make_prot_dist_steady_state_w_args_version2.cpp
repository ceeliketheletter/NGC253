//
//  make_prot_dist_steady_state_new.cpp
//  ICModeling
//
//  Created by Tonia Venters on 12/11/13.
//
//

#include "math.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "TF1.h"
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include <string>

using namespace std;

const double pi=3.14159265358979323846264338327950288;
const double kbeVperK = 8.6173325e-5;
const double elecmass = 5.10999e5; // eV/c^2
const double protonmass = 9.38272013e8; // eV/c^2
const double muonmass = 1.05658367e8; // eV/c^2
const double pionmass = 1.3957018e8; // eV/c^2
const double pi0mass = 1.349766e8; // eV/c^2
const double etamass = 5.47853e8; // eV/c^2
const double elecmass_grams = 9.10938215e-28;
const double protonmass_grams = 1.672621637e-24;
const double cmbtemp = 2.725; // K
const double kBcmbtemp = kbeVperK*cmbtemp; // eV
const double bgphotenergy = 2.7*kbeVperK*cmbtemp; // eV
const double speedoflight = 2.9979246e10;  // cm/s
const double hplank = 4.135667e-15;  // eV*s
const double hbar = hplank/(2.*pi);  // eV*s
const double hc = hplank*speedoflight; // eV*cm
const double thomsoncs = 6.6524e-25;  // cm^2
const double eV2ergs = 1.602e-12;
const double ergs2eV = 1./eV2ergs;
const double pc = 3.0857e18;  // cm
const double kpc = 3.0857e21; // cm
const double mpc = 3.0857e24;  // cm
const double mpc3 = pow(mpc,3.0); // cm^3
const double alphafine = 7.2973525698e-3;
const double year = 525600*60; // s
const double excit_pot_hyd = 14.9913; // eV -- Excitation values from Kamakura et al. 2006, J. Appl. Phys., 100, 064905
const double excit_pot_h2 = 19.2; // eV
const double excit_pot_he = 41.8; // eV
const double atom_mass_hyd = 1.00794; // g mol^-1
const double atom_mass_he = 4.002602; // g mol^-1
const double avo_number = 6.02214129e23; // mol^-1
const double class_elec_radius = 2.1879403267e-13; // cm
const double const_K = 4.*pi*avo_number*class_elec_radius*class_elec_radius*elecmass; // eV mol^-1 cm^2
const double amu = 1.660538921e-24; // g
const double mass_hyd_grams = atom_mass_hyd*amu; // g
const double mass_he_grams = atom_mass_he*amu; // g
const double atom_no_hyd = 1.0;
const double atom_no_he = 2.0;
const double barn = 1.e-24; // cm^2
const double millibarn = (1.e-3)*barn; // cm^2
const double solarmass = 1.9884e33; // g
const double K_pi = 0.17; // related to the fraction of kinetic energy of the proton transferred to gamma rays or leptons -- see Kelner et al. (2006)

// Source parameters

// Please use Gaussian units for constants relating to magnetic field

//const double regBfield = 5.e-5; // Gauss -- range: 50 - 400 microGauss - NGC 253 - Lacki et al. 2012 leptonic
const double regBfield = 4.e-4; // Gauss -- range: 50 - 400 microGauss - NGC 253 - Lacki et al. 2012 hadronic

//const double regBfield = 1.e-5; // Gauss -- 30Dor range: 3-50 (Murphy et al.)
//const double cohlength0 = 1.0; // Mpc -- probably will have to change

const double B_200 = regBfield*1./(200*(1.e-6));

//const double sourcedist = 2.5*mpc; // cm -- NGC 253
const double sourcedist = 3.5*mpc; // cm -- NGC 253 -- Lacki et al. 2012
//const double sourcedistmpc = 2.5; // NGC 253
const double sourcedistmpc = 3.5; // NGC 253 -- Lacki et al. 2012

//const double sourcedist = 54.9*kpc; // cm -- 30Dor
//const double sourcedistmpc = sourcedist*1./mpc; // 30Dor

const double irregion = 150*pc; // radius -- NGC 253
const double irregionmpc = irregion*1./mpc; // NGC 253

//const double irregion = 1000.*pc; // radius -- 30Dor
//const double irregionmpc = irregion*1./mpc; // 30Dor
//const double volirregion = (4./3)*pi*irregion*irregion; // 30Dor

//const double irscaleheight = 70*pc; // NGC 253 Scale height is measured from the midplane to the edge
//const double irscaleheight = 50.0*pc; // NGC 253 -- Lacki et al. 2012 leptonic
const double irscaleheight = 5.0*pc; // NGC 253 -- Lacki et al. 2012 hadronic
const double h_100 = irscaleheight*1./(100.*pc); // NGC 253
const double volirregion = 2.*pi*irregion*irregion*irscaleheight; // NGC 253 -- Disk volume

const double sourcered = 0.0;
//const double wind_speed = 300; // km/s -- NGC 253
const double wind_speed = 1000; // km/s -- NGC 253
//const double v_wind_500 = wind_speed*1./500.; // NGC 253
const double v_wind_300 = wind_speed*1./300.; // NGC 253
//const double nhion = /*600*/1; // cm^-3 -- C&F use 600 cm^-3; Check other groups, as well
//const double nhneut = /*600*/1; // cm^-3
/*const double nh = 1; // cm^-3
 const double gas_no_dens = nh*1./0.9;
 const double nhe = gas_no_dens - nh;*/

//const double gas_mass = 2.e7*solarmass; // g -- range: (2-5)*10^7 M_solar -- NGC 253
const double gas_mass = 3.125e7*solarmass; // g -- range: (2-5)*10^7 M_solar -- NGC 253 -- Lacki et al. 2012
const double gas_density = gas_mass*1./volirregion; // g cm^-3 -- NGC 253

const double f_nh = .9;
const double f_nhe = .1;
const double f_nh_neut = 1.0;
const double f_nh_ion = (1.-f_nh_neut);
const double f_nh_neut_atom = 0;
const double f_nh_neut_mol = (1.-f_nh_neut_atom);

const double nh = gas_density*f_nh*1./mass_hyd_grams; // cm^-3 -- NGC 253
const double nhe = gas_density*f_nhe*1./mass_he_grams; // cm^-3 -- NGC 253
const double n_ism = nh + nhe; // cm^-3 -- NGC 253

//const double n_ism = 2.0; // cm^-3 -- 30Dor
//const double nh = n_ism*f_nh; // cm^-3 -- 30Dor
//const double nhe = n_ism*f_nhe; // cm^-3 -- 30Dor

const double n_250 = n_ism*1./(250);
const double nh_neut = nh*f_nh_neut; // cm^-3
const double nh_ion = nh*f_nh_ion; // cm^-3
const double nhi = nh_neut*f_nh_neut_atom; // cm^-3
const double nh_mol = nh_neut*f_nh_neut_mol*1./2.; // cm^-3
//const double nhi = gas_density*f_nh*f_nh_neut*f_atom_h*1./mass_hyd_grams; // cm^-3
//const double nh_mol = gas_density*f_nh*f_nh_neut*f_mol_h*1./(2.*mass_hyd_grams); // cm^-3

//const double snrate = 0.08/year; // s^-1 -- NGC 253
const double snrate = 0.02/year; // s^-1 -- NGC 253 -- Lacki et al. 2012
//const double snrate = 0.15/year; // s^-1 -- 30Dor

const double snenergy = pow(10.0,51.0)*ergs2eV; // eV
//const double sneff_proton = 0.1;
//const double sneff_proton = 0.2; // Lacki et al. 2012 leptonic
const double sneff_proton = 0.1; // Lacki et al. 2012 hadronic
//const double nh = 1;

// Proton injection spectrum parameters - power law
const double prot_p = 2.1;//2.1; 2.75;  // cf. Chakraborty & Fields 2013 for starburst galaxies -- 2.1; C&F say that range is: 2.0--2.4; Lacki et al. say the range is 2.0--2.6
//const double elec_norm_pl = ;
// 30Dor range: 2.1--2.8 (Murphy et al.)

// Proton injection spectrum parameters - smoothly broken power law
const double prot_p_faint = 1.8;     // cf. Strong et al. 2010
const double prot_p_bright = 2.25;
const double prot_n = 1.0; // strength of the transition
const double prot_break = 4.e9; // eV
const double prot_high_cutoff = 1.e15; // eV
//const double elec_norm_sbpl = ;

// Range of Proton injection spectrum
const double prot_E_min = protonmass; // eV
const double prot_E_max = /*(1.e6)*protonmass*/1.e20; // eV
//const double prot_E_max = (1.e6)*protonmass; // eV
const double log_prot_gmin = log10(prot_E_min*1./protonmass);   //WHAT WHAT WHAT
const double log_prot_gmax = log10(prot_E_max*1./protonmass);
const double zcr = 1.0;

const double pi0_E_max = K_pi*(prot_E_max - protonmass); // eV
const double log_E_pi0_max = log10(pi0_E_max);
const double pion_E_max = pi0_E_max;
const double log_E_pion_max = log10(pion_E_max);

// Range of electron injection spectrum
//const double elec_E_min = 1.5e6; // eV
const double elec_E_min = elecmass;
const double elec_E_max = (1.e6)*elecmass; // eV
const double log_elec_gmin = log10(elec_E_min*1./elecmass);
const double log_elec_gmax = log10(elec_E_max*1./elecmass);

// Spectrum type
// Case 1: Power Law
// Case 2: Power Law with Exponential Cutoff
// Case 3: Smoothly Broken Power Law
// Case 4: Smoothly Broken Power Law with Exponential Cutoff

// const double specflag = 2.0; // Must be 1 or 2 at the moment

// Photon Backgrounds
// Note that all differentials that have to be integrated over energy are to be integrated over log_10(energy) rather than energy. So, if they look a little strange, that would be the reason.

const double isrfmodeldistmpc = 50;  // NGC 253, but you don't need to worry about it here

// Cross sections

double pp_inel_high(double protenergy) // eV
{
    double prot_energy_tev = protenergy*(1.e-12); // TeV
    double ln_prot_energy_tev = log(prot_energy_tev);
    
    //cout << prot_energy_tev << " " << (34.3 + 1.88*ln_prot_energy_tev + 0.25*ln_prot_energy_tev*ln_prot_energy_tev) << endl;
    
    return (34.3 + 1.88*ln_prot_energy_tev + 0.25*ln_prot_energy_tev*ln_prot_energy_tev)*millibarn; //cm^2
}
double pp_inel_low(double protenergy) // eV
{
    double prot_energy_tev = protenergy*(1.e-12); // TeV
    double thr_energy = 1.22e-3; // TeV
    double thr_prot_ratio = thr_energy*1./prot_energy_tev;
    double thr_prot_ratio_sqrd = thr_prot_ratio*thr_prot_ratio;
    double one_min_ratio_fourth = 1.-thr_prot_ratio_sqrd*thr_prot_ratio_sqrd;
    
    //cout << protenergy << " " << pp_inel_high(protenergy)*one_min_ratio_fourth*one_min_ratio_fourth << endl;
    
    return pp_inel_high(protenergy)*one_min_ratio_fourth*one_min_ratio_fourth; // cm^2
}

// Secondary Particle Spectra Per Proton of gamma = gamma_p

// -- Pion decay --
// ---- Emissivities from Kelner, Aharonian, & Bugayov (2006) ----
// ------ E_p > 0.1 TeV ------

double f_gamma_high(double x_gamma, double protenergy) // x_gamma = E_gamma/E_p; [protenergy] = eV
{
    double prot_energy_tev = protenergy*(1.e-12); // TeV
    double ln_p_energy_tev = log(prot_energy_tev);
    double ln_x = log(x_gamma);
    double B_gamma, beta_gamma, k_gamma;
    double x_pow_beta, factor1, factor1_pow_4, factor2;
    
    B_gamma = 1.30 + 0.14*ln_p_energy_tev + 0.011*ln_p_energy_tev*ln_p_energy_tev;
    beta_gamma = 1./(1.79 + 0.11*ln_p_energy_tev + 0.008*ln_p_energy_tev*ln_p_energy_tev);
    k_gamma = 1./(0.801 + 0.049*ln_p_energy_tev + 0.014*ln_p_energy_tev*ln_p_energy_tev);
    
    x_pow_beta = pow(x_gamma,beta_gamma);
    
    factor1 = (1. - x_pow_beta)*1./(1. + k_gamma*x_pow_beta*(1.-x_pow_beta));
    factor1_pow_4 = pow(factor1,4.0);
    factor2 = (1./ln_x) - 4.*beta_gamma*x_pow_beta/(1.-x_pow_beta) - 4.*k_gamma*beta_gamma*x_pow_beta*(1.-2.*x_pow_beta)/(1.+k_gamma*x_pow_beta*(1.-x_pow_beta));
    
    return B_gamma*(ln_x*1./x_gamma)*factor1_pow_4*factor2;
}
double f_e_high(double x_e, double protenergy) // x_e = E_e/E_p; [protenergy] = eV
{
    double prot_energy_tev = protenergy*(1.e-12); // TeV
    double ln_p_energy_tev = log(prot_energy_tev);
    double ln_x = log(x_e);
    double B_e, beta_e, k_e;
    double x_pow_beta, factor1, factor1_pow_3, factor2, factor3;
    
    B_e = 1./(69.5 + 2.65*ln_p_energy_tev + 0.3*ln_p_energy_tev*ln_p_energy_tev);
    beta_e = pow((0.201 + 0.062*ln_p_energy_tev + 0.00042*ln_p_energy_tev*ln_p_energy_tev),-0.25);
    k_e = (0.279 + 0.141*ln_p_energy_tev + 0.0172*ln_p_energy_tev*ln_p_energy_tev)*1./(0.3 + (2.3+ln_p_energy_tev)*(2.3+ln_p_energy_tev));
    
    x_pow_beta = pow(x_e,beta_e);
    
    factor1 = 1. + k_e*ln_x*ln_x;
    factor1_pow_3 = pow(factor1,3.0);
    factor2 = x_e*(1. + 0.3/x_pow_beta);
    factor3 = pow(-1.*ln_x,5.0);
    
    return B_e*(factor1_pow_3*1./factor2)*factor3;
}
double f_pi_high(double x_pi, double protenergy) // x_pion = E_pion/E_p; [protenergy] = eV
{
    double prot_energy_tev = protenergy*(1.e-12); // TeV
    double ln_p_energy_tev = log(prot_energy_tev);
    //double B_pi = 5.58 + 0.78*ln_p_energy_tev + 0.10*ln_p_energy_tev*ln_p_energy_tev;
    double a = 3.67 + 0.83*ln_p_energy_tev + 0.075*ln_p_energy_tev*ln_p_energy_tev;
    double B_pi = a + 0.25;
    //double r = 3.1/pow(B_pi,1.5);
    double r = 2.6/sqrt(a);
    //double alpha = 0.89/(sqrt(B_pi)*(1.-exp(-0.33*B_pi)));
    double alpha = 0.98/sqrt(a);
    double x_alpha = pow(x_pi,alpha);
    double factor1,factor2,factor2n,factor2d,factor2fourth,factor3,factor4;
    
    factor1 = 4.*alpha*B_pi*(x_alpha*1./x_pi);
    factor2n = 1. - x_alpha;
    factor2d = (1. + r*x_alpha*factor2n);
    factor2 = factor2n*1./factor2d;
    factor2fourth = factor2*factor2*factor2*factor2;
    factor3 = (1./factor2n + r*(factor2n-x_alpha)*1./factor2d);
    factor4 = sqrt(1.-(pi0mass*1./(x_pi*protenergy)));
    
    if ((pi0mass*1./(x_pi*protenergy)) > 1.0) factor4 = 0;
    
    //cout << protenergy << " " << factor4 << " " << x_pi << " " << (pi0mass*1./(x_pi*protenergy)) << endl;
    
    return factor1*factor2fourth*factor3*factor4;
}
double f_eta_high(double x_eta, double protenergy) // x_eta = E_eta/E_p; [protenergy] = eV
{
    double ln_x = log(x_eta);
    double expr = (0.55 + 0.028*ln_x)*(1.-etamass*1./(x_eta*protenergy))*f_pi_high(x_eta,protenergy);
    
    return expr;
    //cout << "eta: " << x_eta << " " << expr << endl;
    
    //return (0.55 + 0.028*ln_x)*(1.-etamass*1./(x_eta*protenergy))*f_pi_high(x_eta,protenergy);
}
double e_f_pi_high_diff(double *x_pi, double *par)
{
    double protenergy = par[0];
    double E_pi = x_pi[0]*protenergy;
    double expr = f_pi_high(x_pi[0],protenergy);
    
    //cout << protenergy << " " << x_pi[0] << " " << expr << endl;
    
    return E_pi*expr;
    
    //return E_pi*f_pi_high(x_pi[0],protenergy);
}
double e_f_eta_high_diff(double *x_eta, double *par)
{
    double protenergy = par[0];
    double E_eta = x_eta[0]*protenergy;
    double expr = f_eta_high(x_eta[0],protenergy);
    
    //cout << protenergy << " " << x_eta[0] << " " << expr << endl;
    
    return E_eta*expr;
    
    //return E_eta*f_eta_high(x_eta[0],protenergy);
}

// ------ E_p < 0.1 TeV ------

double q_gamma_p_pi0_low(double pionenergy)
// To be multiplied by the proton spectrum (d^2N/dVdE) -- [q_pion] = s^-1; [E_pi] = eV
// Must be normalized by using the condition of continuity of the spectrum at the point E = 100 GeV
// Note that the normalization is different for gamma rays, electrons, and neutrinos.
{
    double E_p = protonmass + pionenergy*1./K_pi;
    
    //cout << E_p << " " << pp_inel_low(E_p)*(speedoflight*nh*1./K_pi) << endl;
    
    return pp_inel_low(E_p)*(speedoflight*nh*1./K_pi); // pions s^-1 protons^-1
}

// Knock-on

double kn_elec_spectrum_tor(double gamma_prim, double gamma_e) // Currently not in usage; only for testing
{
    //double gamma_prim = pow(10.0,loggamma[0]);
    double gam_sqrd_prim = gamma_prim*gamma_prim;
    //double secenergy = par[0], gamma_e = secenergy*1./elecmass;
    double gamma_e_minus_one = gamma_e - 1.;
    double gamma_e_minus_one_sqrd = gamma_e_minus_one*gamma_e_minus_one;
    double sh = elecmass*1./(atom_mass_hyd*protonmass);
    double sh_sqrd = sh*sh;
    double sh_sqrd_plus_one = sh_sqrd + 1.;
    double she = elecmass*1./(atom_mass_he*protonmass);
    double she_sqrd = she*she;
    double she_sqrd_plus_one = she_sqrd + 1.;
    //double z_tar = 1.0;
    double z_h = 1.0;
    double z_he = 2.0;
    double z_prim = 1.0;
    double phi_h, phi_he;
    double Gamma_1_h = 0.5*sh*(gamma_e-1)+sqrt(1.+(0.5)*(1.+sh*sh)*(gamma_e-1)+(0.25)*sh*sh*(gamma_e-1)*(gamma_e-1));
    double Gamma_1_he = 0.5*she*(gamma_e-1)+sqrt(1.+(0.5)*(1.+she*she)*(gamma_e-1)+(0.25)*she*she*(gamma_e-1)*(gamma_e-1));
    
    //double she = elecmass*1./(atom_mass_he*protonmass);
    double constant_h = 2.0*pi*avo_number*class_elec_radius*class_elec_radius*z_prim*z_prim*z_h*1./(atom_mass_hyd);
    double constant_he = 2.0*pi*avo_number*class_elec_radius*class_elec_radius*z_prim*z_prim*z_he*1./(atom_mass_he);
    
    if (gamma_prim >= Gamma_1_h) phi_h = constant_h*(1./(1.-1./gam_sqrd_prim))*((1./gamma_e_minus_one_sqrd) - ((sh*(gamma_prim + sh_sqrd_plus_one*1./(2.*sh)))*1./(gamma_e_minus_one*gam_sqrd_prim)) + sh_sqrd*1./(2.*gam_sqrd_prim));
    else phi_h = 0;
    
    if (gamma_prim >= Gamma_1_he) phi_he = constant_he*(1./(1.-1./gam_sqrd_prim))*((1./gamma_e_minus_one_sqrd) - ((she*(gamma_prim + she_sqrd_plus_one*1./(2.*she)))*1./(gamma_e_minus_one*gam_sqrd_prim)) + she_sqrd*1./(2.*gam_sqrd_prim));
    else phi_he = 0;
    
    //return 1.75*(0.9*phi_h+4.0*0.1*phi_he);
    
    return phi_h; // cm^2 g^-1
    
    //return constant*(1./(1.-1./gam_sqrd_prim))*((1./gamma_e_minus_one_sqrd) - ((sh*(gamma_prim + sh_sqrd_plus_one*1./(2.*sh)))*1./(gamma_e_minus_one*gam_sqrd_prim)) + sh_sqrd*1./(2.*gam_sqrd_prim)); // cm^2 g^-1
    
    //return constant*(1./(1.-1./gam_sqrd_prim))*(1./((gamma_e - 1.)*(gamma_e - 1.)) - (s*(gamma_prim + (s*s+1)*1./(2*s)))*1./((gamma_e - 1.)*gam_sqrd_prim) )
}
double kn_elec_spectrum_pdg(double gamma_prim, double gamma_e) //dN_{\gamma_p}/dtdGamma_e -- Please check over this function.
{
    //double secenergy = par[0], T_e = secenergy - elecmass;
    double T_e = (gamma_e - 1.)*elecmass;
    double prim_energy = gamma_prim*protonmass;
    double prim_beta_sqrd = 1.0-1.0/(gamma_prim*gamma_prim), prim_beta = sqrt(prim_beta_sqrd);
    double prim_beta_gam = prim_beta*gamma_prim;
    double prim_beta_gam_sqrd = prim_beta_gam*prim_beta_gam;
    double z_prim = 1.0;
    double thickness_hi = nhi*mass_hyd_grams*prim_beta*speedoflight; //g cm^-2 s^-1
    double thickness_h2 = nh_mol*(2.*mass_hyd_grams)*prim_beta*speedoflight; // g cm^-2 s^-1
    double thickness_he = nhe*mass_he_grams*prim_beta*speedoflight; //g cm^-2 s^-1
    double constant_hi = elecmass*thickness_hi*(0.5)*const_K*z_prim*z_prim*(atom_no_hyd*1.0/atom_mass_hyd)*(1./prim_beta_sqrd); // eV^2 s^-1
    double constant_h2 = elecmass*thickness_h2*(0.5)*const_K*z_prim*z_prim*(atom_no_hyd*1.0/atom_mass_hyd)*(1./prim_beta_sqrd); // eV^2 s^-1;
    double constant_he = elecmass*thickness_he*(0.5)*const_K*z_prim*z_prim*(atom_no_he*1.0/atom_mass_he)*(1./prim_beta_sqrd); // eV^2 s^-1
    double T_max = 2.*elecmass*(prim_beta_gam_sqrd*1./(1.+2.*gamma_prim*(elecmass*1./protonmass)+(elecmass*1./protonmass)*(elecmass*1./protonmass))); // if using higher Z primaries, must replace protonmass with the mass of the primary
    double f_of_t_hi, f_of_t_h2, f_of_t_he;
    
    //const_K = 4.*pi*avo_number*class_elec_radius*class_elec_radius*elecmass; // eV mol^-1 cm^2
    
    if ((T_e > excit_pot_hyd)&&(T_e <= T_max))
        f_of_t_hi = 1. - prim_beta_sqrd*(T_e*1./T_max) + 0.5*(T_e*1./prim_energy)*(T_e*1./prim_energy); // spin 1/2 formula -- other formulae are given in the notes and in High Energy Particles by Bruno Rossi (1952) or the Particle Data Group Review
    else f_of_t_hi = 0;
    if ((T_e > excit_pot_h2)&&(T_e <= T_max))
        f_of_t_h2 = 1. - prim_beta_sqrd*(T_e*1./T_max) + 0.5*(T_e*1./prim_energy)*(T_e*1./prim_energy);
    else f_of_t_h2 = 0;
    if ((T_e > excit_pot_he)&&(T_e <= T_max))
        f_of_t_he = 1. - prim_beta_sqrd*(T_e*1./T_max) + 0.5*(T_e*1./prim_energy)*(T_e*1./prim_energy);
    else f_of_t_he = 0;
    
    //cout << T_e << " " << T_max << " " << f_of_t << endl;
    
    //return constant*f_of_t*1./(T_e*T_e); // s^-1
    return constant_hi*f_of_t_hi*1./(T_e*T_e) + constant_h2*f_of_t_h2*1./(T_e*T_e) + constant_he*f_of_t_he*1./(T_e*T_e); // s^-1
}

// Primary Particle Continuous Losses

double bE_pp_high(double protenergy)
{
    double epiinteg,eetainteg;
    
    int np = 25;
    double *x = new double[np];
    double *w = new double[np];
    
    TF1 e_f_pi_h_diff("e_f_pi_h_diff",e_f_pi_high_diff,-0.5,1.5,1);
    e_f_pi_h_diff.SetParameter(0,protenergy);
    e_f_pi_h_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    epiinteg = e_f_pi_h_diff.IntegralFast(np,x,w,(1.0e-3),1.0);
    
    delete [] x;
    delete [] w;
    
    double *x1 = new double[np];
    double *w1 = new double[np];
    
    TF1 e_f_eta_h_diff("e_f_eta_h_diff",e_f_eta_high_diff,-0.5,1.5,1);
    e_f_eta_h_diff.SetParameter(0,protenergy);
    e_f_eta_h_diff.CalcGaussLegendreSamplingPoints(np,x1,w1,1.e-15);
    eetainteg = e_f_eta_h_diff.IntegralFast(np,x1,w1,(1.0e-3),1.0);
    
    delete [] x1;
    delete [] w1;
    
    return speedoflight*nh*pp_inel_high(protenergy)*(epiinteg+eetainteg); // eV s^-1
}
double bE_pp_low(double protenergy)
{
    double prot_ke = protenergy - protonmass;
    double pionenergy = K_pi*prot_ke;
    return speedoflight*nh*pp_inel_low(protenergy)*pionenergy;
}
double bE_pp(double loggamma)
{
    double protenergy = pow(10.0,loggamma)*protonmass;
    double protenergy_gev = protenergy*(1.e-9);
    double normalization, be_high, be_low;
    
    normalization = bE_pp_high(1.e11)*1./bE_pp_low(1.e11);
    
    be_low = normalization*bE_pp_low(protenergy);
    be_high = bE_pp_high(protenergy);
    
    if (protenergy_gev >= 1.22){
        if (protenergy <= (1.e11)) return be_low;
        else return be_high;
    }
    else return 0.0;
}
double bE_pp_tor(double loggamma) // Currently not in use; for testing purposes only
{
    double gamma = pow(10.0,loggamma);
    double protenergy_gev = gamma*protonmass*(1.e-9); // GeV
    double prot_ke_gev = (gamma - 1.)*protonmass*(1.e-9); // GeV
    double sigma_pp_inel = 3.e-26; // cm
    
    if (protenergy_gev < 100) sigma_pp_inel = pp_inel_low(protenergy_gev*(1.e9));
    else sigma_pp_inel = pp_inel_high(protenergy_gev*(1.e9));
    
    //cout << sigma_pp_inel << endl;
    
    //if (protenergy_gev >= 1.22) return (5.85e-16)*nh*prot_ke_gev*(1.e9); // eV s^-1
    //if (protenergy_gev >= 1.22) return 0.65*speedoflight*sigma_pp_inel*nh*prot_ke_gev*(1.e9); // eV s^-1
    if (protenergy_gev >= 1.22) return 0.17*speedoflight*sigma_pp_inel*nh*prot_ke_gev*(1.e9); // eV s^-1
    else return 0.0;
}
double bE_prot_synch(double loggamma, double bfield)
{
    
    double gamma = pow(10.0,loggamma);
	double gammasqrd = gamma*gamma;
    double zprim = 1.0;
    double zprim_fourth = pow(zprim,4.0);
    double constant = (4.0/3)*zprim_fourth*thomsoncs*speedoflight*(elecmass*1./protonmass)*(elecmass*1./protonmass); // if using higher Z primaries, must replace protonmass with the mass of the primary
    //double Bfield = B0*(1.+redshift);
	double bfieldsqrd = bfield*bfield;
	double magenergydensity = (bfieldsqrd/(8.*pi))*ergs2eV; // eV cm^-3
    
    //cout << magenergydensity << endl;
    
    return constant*(gammasqrd-1.)*magenergydensity; //(gamma^2-1) = gamma^2*beta^2 // eV s^-1
}
double bE_ion_tor(double loggamma) // Off by a factor of two?
{
    double gamma = pow(10.0,loggamma);
    double beta_sqrd = 1.-1./(gamma*gamma);
    double beta = sqrt(beta_sqrd);
    //double pr_energy = gamma*protonmass;
    
    return (1.83e-17)*(1.e9)*nh*(1./beta)*(10.9 + 2.*log(gamma) + log(beta_sqrd) - beta_sqrd); // eV s^-1
}
double bE_ion_pdg(double loggamma)
{
    double gamma = pow(10.0,loggamma);
    double beta_sqrd = 1.-1./(gamma*gamma);
    double beta = sqrt(beta_sqrd);
    double beta_gamma = beta*gamma;
    double beta_gamma_sqrd = beta_gamma*beta_gamma;
    double z_prim = 1.0;
    double thickness_hi = nhi*mass_hyd_grams*beta*speedoflight; // g cm^-2 s^-1
    double thickness_h2 = nh_mol*(2.*mass_hyd_grams)*beta*speedoflight; // g cm^-2 s^-1
    double thickness_he = nhe*mass_he_grams*beta*speedoflight; // g cm^-2 s^-1
    double constant_hi = thickness_hi*const_K*z_prim*z_prim*(atom_no_hyd*1.0/atom_mass_hyd)*(1./beta_sqrd); // eV s^-1
    double constant_h2 = thickness_h2*const_K*z_prim*z_prim*(atom_no_hyd*1.0/atom_mass_hyd)*(1./beta_sqrd); // eV s^-1
    double constant_he = thickness_he*const_K*z_prim*z_prim*(atom_no_he*1.0/atom_mass_he)*(1./beta_sqrd); // eV s^-1
    double T_max = 2.*elecmass*(beta_gamma_sqrd*1./(1.+2.*gamma*(elecmass*1./protonmass)+(elecmass*1./protonmass)*(elecmass*1./protonmass)));
    double f_of_t_hi = (0.5)*log(2.*elecmass*beta_gamma_sqrd*T_max/(excit_pot_hyd*excit_pot_hyd)) - beta_sqrd;
    double f_of_t_h2 = (0.5)*log(2.*elecmass*beta_gamma_sqrd*T_max/(excit_pot_h2*excit_pot_h2)) - beta_sqrd;
    double f_of_t_he = (0.5)*log(2.*elecmass*beta_gamma_sqrd*T_max/(excit_pot_he*excit_pot_he)) - beta_sqrd;
    
    //return constant*f_of_t; // eV s^-1
    return constant_hi*f_of_t_hi + constant_h2*f_of_t_h2 + constant_he*f_of_t_he; // eV s^-1
}
double bE_total(double loggamma, double bfield)
{
    double gamma = pow(10.0, loggamma);
    //double beta = sqrt(1.-1./(gamma*gamma));
    
    //cout << gamma*protonmass << " " << bE_pp(loggamma) << " " << bE_prot_synch(loggamma,regBfield) << " " << bE_ion_pdg(loggamma) << endl;
    
    return (bE_pp(loggamma)+bE_prot_synch(loggamma,bfield)+bE_ion_pdg(loggamma)); // eV s^-1
}

// -- Primary Injection Spectra -- dN/dEdVdt

double powerlaw(double energy, double power)
{
    //double power = par[0];
    return pow(energy,-1.0*power);
}
double cutpowerlaw(double energy, double *par) // No. parameters: 2
{
    double power = par[0];
    double ecut = par[1];
    
    //cout << power << " " << ecut << endl;
    
    return pow(energy,-1.0*power)*exp(-1.*energy/ecut);
}
double smoothbpl(double energy, double *par) // No. parameters: 4
{
    double faint_slope = par[0];
    double bright_slope = par[1];
    double n = par[2];
    double breakenergy = par[3];
    
    double factor1 = pow((energy*1./breakenergy),faint_slope*n);
    double factor2 = pow((energy*1./breakenergy),bright_slope*n);
    
    return pow((factor1+factor2),-1.0/n);
}
double cutsmoothbpl(double energy, double *par) // No. parameters: 5
{
    double faint_slope = par[0];
    double bright_slope = par[1];
    double n = par[2];
    double breakenergy = par[3];
    double ecut = par[4];
    
    double factor1 = pow((energy*1./breakenergy),faint_slope*n);
    double factor2 = pow((energy*1./breakenergy),bright_slope*n);
    
    return pow((factor1+factor2),-1.0/n)*exp(-1.*energy/ecut);
}

// -- Steady state solutions
// -- -- Infinite escape timescale

double QE_no_escape_diff(double *loggamma, double *par) // No. parameters: 1 // Parameters checked
{
    double qe;
    double energy = pow(10.0,loggamma[0])*protonmass;
    int flag = (int)(par[0]+0.5);
    double bpl_par[5];
    
    //cout << flag << endl;
    
    /*switch (flag) {
        case 1: // Power Law
            qe = powerlaw(energy,prot_p)*energy*log(10.0);
            //cout << "Hey!\n";
            break;
        case 2: // Cutoff Power Law
            //double bpl_par[2];
            //bpl_par.push_back(elec_p);
            //bpl_par.push_back(elec_high_cutoff);
            
            bpl_par[0] = prot_p;
            bpl_par[1] = prot_high_cutoff;
            
            qe = cutpowerlaw(energy,bpl_par)*energy*log(10.0);
            break;
        case 3: // Smoothly Broken Power Law
            //double bpl_par[4];
            
            bpl_par[0] = prot_p_faint;
            bpl_par[1] = prot_p_bright;
            bpl_par[2] = prot_n;
            bpl_par[3] = prot_break;
            
            qe = smoothbpl(energy,bpl_par)*energy*log(10.0);
            break;
        case 4: // Cutoff Smoothly Broken Power Law
            //double bpl_par[5];
            
            bpl_par[0] = prot_p_faint;
            bpl_par[1] = prot_p_bright;
            bpl_par[2] = prot_n;
            bpl_par[3] = prot_break;
            bpl_par[4] = prot_high_cutoff;
            
            qe = cutsmoothbpl(energy,bpl_par)*energy*log(10.0);
            break;
        default:
            qe = powerlaw(energy,prot_p)*energy*log(10.0);
            break;
    }*/
    
    //Spectral parameters -- inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut
    switch (flag) {
        case 1: // Power Law
            qe = powerlaw(energy,par[1])*energy*log(10.0);
            //cout << "Hey!\n";
            break;
        case 2: // Cutoff Power Law
            bpl_par[0] = par[1];
            bpl_par[1] = par[5];
            
            //cout << cutpowerlaw(energy,bpl_par) << endl;
            
            qe = cutpowerlaw(energy,bpl_par)*energy*log(10.0);
            break;
        case 3: // Smoothly Broken Power Law
            bpl_par[0] = par[1];
            bpl_par[1] = par[2];
            bpl_par[2] = par[3];
            bpl_par[3] = par[4];
            
            qe = smoothbpl(energy,bpl_par)*energy*log(10.0);
            break;
        case 4: // Cutoff Smoothly Broken Power Law
            bpl_par[0] = par[1];
            bpl_par[1] = par[2];
            bpl_par[2] = par[3];
            bpl_par[3] = par[4];
            bpl_par[4] = par[5];
            
            qe = cutsmoothbpl(energy,bpl_par)*energy*log(10.0);
            break;
        default:
            qe = powerlaw(energy,par[1])*energy*log(10.0);
            break;
    }
    
    //cout << energy << " " << qe << endl;
    
    return qe; // Q(E) = dN/d(log(gamma))dVdt
}
double QE_kin_diff(double *loggamma, double *par)
{
    double gamma = pow(10.0,loggamma[0]);
    double kinenergy = (gamma - 1)*protonmass;
    
    //Parameters -- specpar,inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut
    
    //cout << kinenergy << " " << QE_no_escape_diff(loggamma,par)*kinenergy << endl;
    
    return QE_no_escape_diff(loggamma,par)*kinenergy; // eV cm^-3 s^-1
}
double NE_no_escape(double loggamma, double bfield, int specflag)
{
    double energy = pow(10.0,loggamma)*protonmass; // eV
    double qeinteg;
    //double specflag = 1.0; // Must be 1 or 2 at the moment
    
    //cout << "Here!\n";
    
    int np = 25;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 qe_simple_diff("qe_simple_diff",QE_no_escape_diff,(log_prot_gmin-0.5),(log_prot_gmax+0.5),1);
    qe_simple_diff.SetParameter(0,specflag);
    qe_simple_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qeinteg = qe_simple_diff.IntegralFast(np,x,w,loggamma,log_prot_gmax); // dN/dVdt
    
    delete [] x;
	delete [] w;
    
    //cout << qeinteg << " " << bE_total(loggamma,regBfield) << endl;
    
    //if (qeinteg < 0) cout << loggamma << " " << log_prot_gmax << endl;
    
    //cout << qeinteg << " " << bE_total(loggamma,regBfield,secenergybin) << endl;
    
    //return qeinteg;
    
    return qeinteg*(1./bE_total(loggamma,bfield))*energy*log(10.0); // dN/d(log(gamma))dV
    
    //return qeinteg*(1./bE_total(loggamma,regBfield))*energy*log(10.0); // dN/d(log(gamma))dV
    //return qeinteg*(1./bE_total(loggamma,regBfield)); // dN/dEdV
}

// -- Finite Escape Timescale

double diff_timescale(double energy) // [energy] = eV
{
    double gamma = energy*1./protonmass;
    double beta_sqrd = 1.0-1.0/(gamma*gamma), beta = sqrt(beta_sqrd);
    
    // NGC 253
    
    //return 2.6e7*year/sqrt(energy*1./(3.e9));  // Milky Way
    //return 1.0e7*beta*year/sqrt(energy*1./(1.e9));  // Paglione & Abrahams 2012
    //return 3.5e6*year*h_100*(1./B_200)*sqrt(n_250); // Lacki et al. 2011
    //return 1.0e7*year/sqrt(energy*1./(3.e9)); // Lacki et al. 2012 leptonic
    return 1.0e6*year/sqrt(energy*1./(3.e9)); // Lacki et al. 2012 hadronic
    
    // 30Dor
    
    //return 3.0e6*year; // Murphy et al. 2012
}
double adv_timescale(double energy) // [energy] = eV
{
    //return 2.0e5*year/v_wind_500; // Lacki et al. 2011
    //return 1.0e6*year; // Paglione & Abrahams 2012
    return 3.0e5*year*h_100*1./v_wind_300;
}
double conv_timescale(double energy)
{
    double tau_diff = diff_timescale(energy);
    double tau_adv = adv_timescale(energy);
    double inv_tau_conv = 1./tau_diff + 1./tau_adv;
    
    return 1./inv_tau_conv;
}
double integ_factor(double *loggamma, double *par)
{
    double gamma = pow(10.0,loggamma[0]);
    double primenergy = gamma*protonmass;
    double factor = primenergy*log(10.0);
    double t_diff = diff_timescale(primenergy);
    double t_conv = conv_timescale(primenergy);
    //double secenergybin = par[0];
    
    //cout << primenergy << " " << bE_total(loggamma[0],regBfield,par[0]) << endl;
    
    //cout << primenergy << " " << t_diff << " " << bE_total(loggamma[0],regBfield) << endl;
    
    return factor*1./(t_diff*bE_total(loggamma[0],regBfield)); // Lacki et al. 2012 hadronic
    //return factor*1./(t_conv*bE_total(loggamma[0],regBfield)); // Lacki et al. 2012 leptonic
}
double QE_mod(double *loggamma, double *par) // par[0] = log10(secprotonenergy)
{
    double gamma = pow(10.0,loggamma[0]);
    double primenergy = gamma*protonmass;
    //double secprotenergy = par[0];
    //double secprotgamma = secprotenergy*1./protonmass;
    //double logsecprgamma = log10(secprotgamma);
    double logsecprgamma = par[0];
    //double plindex = par[1];
    //double logEcut = par[2];
    //double factor = primenergy*log(10.0);
    double integ, exp_factor;
    double qe_par[6];
    
    //specpar,inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut
    
    //qe_par[0] = specflag;
    qe_par[0] = par[1];
    qe_par[1] = par[2];
    qe_par[2] = par[3];
    qe_par[3] = par[4];
    qe_par[4] = par[5];
    qe_par[5] = par[6];
    
    int np = 25;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 integ_fac("integ_fac",integ_factor,(log_prot_gmin-0.5),(log_prot_gmax+0.5),1);
    integ_fac.SetParameter(0,0);
    integ_fac.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    integ = integ_fac.IntegralFast(np,x,w,logsecprgamma,loggamma[0]);
    
    delete [] x;
	delete [] w;
    
    exp_factor = exp(-1.*integ);
    
    //if (exp_factor > 1) cout << integ << endl;
    
    //cout << integ << " " << exp_factor << endl;
    
    //exp_factor = 1.0;
    
    return QE_no_escape_diff(loggamma,qe_par)*exp_factor;
}
double NE_finite_time(double *loggamma, double *par) // par[0] = normalization from main()
//double NE_finite_time(double loggamma)
{
    double gamma = pow(10.0,loggamma[0]);
    double energy = gamma*protonmass;
    //double qeinteg, integ_intfac, expfac, factor = gamma*protonmass*log(10.0);
    double qeinteg;
    double norm = par[0];
    double magfield = par[1];
    //double plindex = par[2];
    //double logEcut = par[3];
    double ne_ss;
    //double specflag = 1.0; // Must be 1 or 2 at the moment
    //normalization,galBfield,specpar,inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut
    
    double funcpars[7];
    
    funcpars[0] = loggamma[0];
    funcpars[1] = par[2];
    funcpars[2] = par[3];
    funcpars[3] = par[4];
    funcpars[4] = par[5];
    funcpars[5] = par[6];
    funcpars[6] = par[7];
    
    int np = 25;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 qe_mod_diff("qe_mod_diff",QE_mod,(log_prot_gmin-0.5),(log_prot_gmax+0.5),7);
    //qe_mod_diff.SetParameters(loggamma[0],plindex,logEcut);
    qe_mod_diff.SetParameters(funcpars);
    qe_mod_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qeinteg = qe_mod_diff.IntegralFast(np,x,w,loggamma[0],log_prot_gmax); // dN/dVdt
    
    delete [] x;
	delete [] w;
    
    //ne_ss = norm*qeinteg*(1./bE_total(loggamma[0],regBfield))*energy*log(10.0); // dN/dVd(log(gamma))
    ne_ss = norm*qeinteg*(1./bE_total(loggamma[0],magfield))*energy*log(10.0); // dN/dVd(log(gamma))
    
//    cout << loggamma[0] << " " << ne_ss << endl;
    
    return ne_ss;
    
    //return qeinteg*(1./bE_total(loggamma,regBfield))*energy*log(10.0); // dN/dVd(log(gamma))
    //return norm*qeinteg*(1./bE_total(loggamma[0],regBfield))*energy*log(10.0); // dN/dVd(log(gamma))
    //return qeinteg*(1./bE_total(loggamma[0],regBfield))*energy*log(10.0); // dN/dVd(log(gamma))
}



int main(){
//int main(int argc, char **argv){     
   
    
     
       
  std::string prefix = "finding_pf_"; string suffix = ".txt"; 
  int i, j, k;
    
  j = 0;
  k =0;
        

    
   
  double galBfield = 4.e-4;   //define non-changing variables outside for loop //remove if calculating galBfield
//  double galBfield;                                                        //add if calculating galBfield
  double specpar = 2.0;
  double inputplleindex;
//  double inputplleindex = 2.1;
  double inputbplhiindex = 0;
  double inputbpln = 0;
  double inputbplEbreak = pow(10.0,0);
  double exponent;
  double ne_sspars[8];
//  double inputlogEcut;
  double inputlogEcut = pow(10.0, 15);    //add if not activating 3rd for loop
  
  double logNnotprime, t, b, y;
  
  
 

  
 
  
  
  
  
  int specflag = (int)(specpar+0.5);
    
        double logbinsize = .1;
        double loggamma = log_prot_gmin + logbinsize;
        double loggamma_0 = 0.1;
        double ne_ss, normalization = 1.0, qekininteg;
        unsigned int num = 0;
        char rootfilename[250];
        char funcname[250];
  
   ne_sspars[0] = normalization;
   ne_sspars[2] = specpar;
   ne_sspars[4] = inputbplhiindex;
   ne_sspars[5] = inputbpln;
   ne_sspars[6] = inputbplEbreak;
   ne_sspars[7] = inputlogEcut;  //remove if activating 3rd for loop
   ne_sspars[1] = galBfield; // remove if activating 2nd loop
 

    
    
  for (i=0, inputplleindex=1.8; inputplleindex<=2.7; inputplleindex+=0.025, ++i)
  {  
    ne_sspars[3] = inputplleindex;
/*    for (j=0, galBfield=5.e-5; galBfield<=4.e-4; galBfield+=0.00007, ++j)    //current range is 5.e-5 -- 4.e-4 Gauss  //add if calculating galBfield
    {
      ne_sspars[1] = galBfield; */
/*      for (k=0, exponent=15; exponent<=20; ++exponent, ++k)    //remove if calculating galBfield
        {
        inputlogEcut = pow(10.0, exponent);     //remove if calculating galBfield
        ne_sspars[7] = inputlogEcut;
*/        

        

        std::string file = prefix + std::to_string (i) + std::to_string (j) + std::to_string (k) + suffix; //concatenates the ith round to the file name
        



        //TGraph protdist;
        //ofstream fne("/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_exp_cutoff.dat");
    
        // Parameters to be input on the command line
        // galBfield [currently set to 4.e-4; current range is 5.e-5 -- 4.e-4 Gauss]
    
        // Specpar: [1.0] Power law; [2.0] Power Law with Exponential Cutoff; [3.0] Smoothly Broken Power Law; [4.0] Smoothly Broken Power Law with Exponential Cutoff
        // -- For [1.0], use inputplleindex [currently set to 2.1; has to be > 2.0]
        // -- For [2.0], use inputplleindex [currently set to 2.1] and inputlogEcut [currently set to 15.0; at 20.0, this is pretty much a power law]
        // -- Ignore for now: For [3.0], use inputplleindex [currently set to 2.1], inputbplhiindex [--], inputbpln [usually I take 1.0], inputbplEbreak [--]
        // -- Ignore for now: For [4.0], use them all
    
        // To run, type: |> ./<name of executable>.out [galBfield value] [specpar] [inputplleindex] [inputbplhiindex] [inputbpln] [inputbplEbreak] [inputlogEcut]
    
        // To run the current model, type: |> ./<name of executable>.out 400 2.0 2.1 0.0 0.0 0.0 15.0
    
        //double galBfield = atof(argv[1]), inputplindex = atof(argv[2]), inputlogEbreak = atof(argv[3]);
//        double galBfield = (atof(argv[1])*(1.e-6)), specpar = atof(argv[2]), inputplleindex = atof(argv[3]), inputbplhiindex = atof(argv[4]), inputbpln = atof(argv[5]), inputbplEbreak = pow(10.0,atof(argv[6])), inputlogEcut = pow(10.0,atof(argv[7]));
    
        switch (specflag) 
        {
            case 1: // Power Law
                sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_plindex_%.1f.root",(galBfield*1.e6),inputplleindex);
                break;
            case 2: // Cutoff Power Law
                sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_plindex_%.1f_Ecut_%.1f.root",(galBfield*1.e6),inputplleindex,inputlogEcut);
                break;
            case 3: // Smoothly Broken Power Law
                sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_leindex_%.1f_hiindex_%.1f_n_%.1f_Eb_%.1f.root",(galBfield*1.e6),inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak);
                break;
            case 4: // Cutoff Smoothly Broken Power Law
                sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_leindex_%.1f_hiindex_%.1f_n_%.1f_Eb_%.1f_Ecut_%.1f.root",(galBfield*1.e6),inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut);
                break;
            default:
                sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_plindex_%.1f.root",(galBfield*1.e6),inputplleindex);
                break;
        }
    
    
    
    
    
        ofstream myfile;
        myfile.open (file);
    
        myfile << "# " << "Output file of: 'make_prot_dist_steady_state_w_args_updated.cpp' " << "\t" << "Round " << i << j << k << endl;
        myfile << "# " << "for galBfield = " << galBfield << "," << " specpar = " << specpar << "," <<" inputplleindex = " << inputplleindex << "," << " and inputlogEcut = " << inputlogEcut << endl;
//        myfile << "Given : " << "[" << argv[1] << "  " << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5] << "  " << argv[6] << "  " << argv[7] << "]" <<endl;
        myfile << '#' << '\n';
        myfile << "# " << "loggamma" << "\t" << "t" <<"\t" << "\t" << "y" << "\t" <<"\t" << "b" << "\t" << "\t" << "inputplleindex"<< endl;;
    
    

        
                    
        for (loggamma=.1; loggamma<=3; loggamma+=logbinsize)
            {
                ne_ss= NE_finite_time(& loggamma, ne_sspars);
    
    
        
                t = loggamma - 0.1;
                logNnotprime = log10(NE_finite_time (&loggamma_0, ne_sspars));
                b = logNnotprime + t;
                y = log10(ne_ss); 
                
                
   
                
                myfile << loggamma << "\t" << "\t" << t << "\t" << "\t" << y << "\t" << "\t" << b << "\t" << "\t" << inputplleindex << endl; //to save to a file, use myfile instead of cout   
                
                                   
                                                        
            }
        myfile.close();
            
            
            
            
//        }    
//    }  
  }  
    
    
    
    //sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_dist_%.1f_scaleheight_%.1f_wind_%.1f_gasmass_%.2f.root",(regBfield*1.e6),sourcedistmpc,(irscaleheight*1./pc),wind_speed,(gas_mass*1.e-6/solarmass));
    //sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_plindex_%.1f_Ecut_%.1f.root",(galBfield*1.e6),inputplindex,inputlogEbreak);
//    sprintf(funcname,"ne_prot_log_prot_gmax_%.1f_eff_%.1f",log_prot_gmax,sneff_proton);
    
    //sprintf(rootfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%d_%d_sourcedist_%d_%d_region_%d_%d_scaleheight_%d_%d_wind_%d_%d_gasmass_%d_%d_f_He_%d_%d_SNrate_%d_%d.root",int(trunc(regBfield*1.e6)),int(((regBfield*1.e6)-trunc(regBfield*1.e6))*10.0),int(trunc(sourcedistmpc)),int((sourcedistmpc-trunc(sourcedistmpc))*10.0),int(trunc(2.0*irregion/pc)),int(((2.0*irregion/pc)-trunc(2.0*irregion/pc))*10.0),int(trunc(irscaleheight*1./pc)),int((irscaleheight*1./pc-trunc(irscaleheight*1./pc))*10.0),int(trunc(wind_speed)),int((wind_speed-trunc(wind_speed))*10.0),int(trunc(gas_mass*1.e-6/solarmass)),int((gas_mass*1.e-6/solarmass-trunc(gas_mass*1.e-6/solarmass))*100.0),int(trunc(f_nhe)),int((f_nhe-trunc(f_nhe))*10.0),int(trunc(snrate*year)),int((snrate*year-trunc(snrate*year))*100.0));
    //sprintf(graphname,"ne_prot_pl_exp_index_%d_%d_Ecut_%d_%d_log_prot_gmax_%d_%d_eff_%d_%d",int(trunc(prot_p)),int((prot_p-trunc(prot_p))*10.0),int(trunc(log10(prot_high_cutoff))),int((log10(prot_high_cutoff)-trunc(log10(prot_high_cutoff)))*10.0),int(trunc(log_prot_gmax)),int((log_prot_gmax-trunc(log_prot_gmax))*10.0),int(trunc(sneff_proton)),int((sneff_proton-trunc(sneff_proton))*10.0));
    
//    TFile protdistrootfile(rootfilename,"recreate");
    
    // Find normalization for proton injection spectrum
    
    /*int np = 25;
	double *x = new double[np];
	double *w = new double[np];
    
    //TF1 *qe_norm_diff = new TF1("qe_kin_diff",QE_kin_diff,0.0,(log_elec_gmax+0.5),1);
    //TF1 *qe_norm_diff = new TF1("qe_kin_diff",QE_kin_diff,log_prot_gmin-.5,log_prot_gmax+.5,1);
    TF1 qe_norm_diff("qe_kin_diff",QE_kin_diff,log_prot_gmin-.5,6.5,6);
    qe_norm_diff.SetParameters(specpar,inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut);
    qe_norm_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qekininteg = qe_norm_diff.IntegralFast(np,x,w,log_prot_gmin,6.0); // dT/dVdt
    //qekininteg = qe_norm_diff->IntegralFast(np,x,w,log_prot_gmin,log_prot_gmax); // dT/dVdt
    //qekininteg = qe_norm_diff->IntegralFast(np,x,w,0.0,100.0); // dT/dVdt
    
    delete [] x;
	delete [] w;
    
    normalization = sneff_proton*snenergy*snrate*1./(volirregion*qekininteg); // dimensionless
    
    cout << normalization << endl;*/
    
//    TF1 protdist(funcname,NE_finite_time,log_prot_gmin+logbinsize,log_prot_gmax,8);
//    protdist.SetParameters(normalization,galBfield,specpar,inputplleindex,inputbplhiindex,inputbpln,inputbplEbreak,inputlogEcut);
    
    //cout << normalization << endl;
    
    /*while (loggamma <= log_prot_gmax) {
        ne_ss = normalization*NE_finite_time(loggamma);
        //ne_ss = normalization*NE_no_escape(loggamma);
        
        cout << loggamma << " " << ne_ss << endl;
        protdist.SetPoint(num,loggamma,ne_ss);
        //fne << loggamma << " " << ne_ss << endl;
        
        loggamma += logbinsize;
        num++;
    }*/
    
    //protdist.SetName(funcname);
//    protdist.Write();
    
//    protdistrootfile.Close();
    
    //fne.close();
//    return 0;
}