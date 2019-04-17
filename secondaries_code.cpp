//
//  make_prot_qe_dist_steady_state_new.cpp
//  
//
//  Created by Tonia Venters on 12/11/13.
//
//

#include "math.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include <vector>
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

//const double regBfield = 5.e-5; // Gauss -- range: 50 - 400 microGauss - NGC 253
double regBfield = 4.e-4; // Gauss -- range: 50 - 400 microGauss - NGC 253

//const double regBfield = 1.e-5; // Gauss -- 30Dor range: 3-50 (Murphy et al.)
//const double cohlength0 = 1.0; // Mpc -- probably will have to change

//const double B_200 = regBfield*1./(200*(1.e-6));

//const double sourcedist = 2.5*mpc; // cm -- NGC 253
//const double sourcedistmpc = 2.5; // NGC 253
double sourcedist = 3.5*mpc; // cm -- NGC 253 -- Lacki et al. 2012
double sourcedistmpc = 3.5; // NGC 253

//const double sourcedist = 54.9*kpc; // cm -- 30Dor
//const double sourcedistmpc = sourcedist*1./mpc; // 30Dor

double irregion = 150*pc; // radius -- NGC 253
double irregionmpc = irregion*1./mpc; // NGC 253

//const double irregion = 1000.*pc; // radius -- 30Dor
//const double irregionmpc = irregion*1./mpc; // 30Dor
//const double volirregion = (4./3)*pi*irregion*irregion; // 30Dor

//const double irscaleheight = 70*pc; // NGC 253 Scale height is measured from the midplane to the edge
//const double irscaleheight = 50.0*pc; // NGC 253 -- Lacki et al. 2012 leptonic
double irscaleheight = 5.0*pc; // NGC 253 -- Lacki et al. 2012 hadronic
double h_100 = irscaleheight*1./(100.*pc); // NGC 253
double volirregion = 2.*pi*irregion*irregion*irscaleheight; // NGC 253 -- Disk volume

const double sourcered = 0.0;
//const double wind_speed = 300; // km/s -- NGC 253
//const double v_wind_500 = wind_speed*1./500.; // NGC 253
double wind_speed = 1000; // km/s -- NGC 253 -- Lacki et al. 2012
double v_wind_300 = wind_speed*1./300.; // NGC 253
//const double nhion = /*600*/1; // cm^-3 -- C&F use 600 cm^-3; Check other groups, as well
//const double nhneut = /*600*/1; // cm^-3
/*const double nh = 1; // cm^-3
 const double gas_no_dens = nh*1./0.9;
 const double nhe = gas_no_dens - nh;*/

//const double gas_mass = 2.e7*solarmass; // g -- range: (2-5)*10^7 M_solar -- NGC 253
double gas_mass = 3.125e7*solarmass; // g -- range: (2-5)*10^7 M_solar -- NGC 253 -- Lacki et al. 2012
double gas_density = gas_mass*1./volirregion; // g cm^-3 -- NGC 253

double f_nh = .9;
double f_nhe = .1;
double f_nh_neut = 1.0;
double f_nh_ion = (1.-f_nh_neut);
double f_nh_neut_atom = 0;
double f_nh_neut_mol = (1.-f_nh_neut_atom);

double nh = gas_density*f_nh*1./mass_hyd_grams; // cm^-3 -- NGC 253
double nhe = gas_density*f_nhe*1./mass_he_grams; // cm^-3 -- NGC 253
double n_ism = nh + nhe; // cm^-3 -- NGC 253

//const double n_ism = 2.0; // cm^-3 -- 30Dor
//const double nh = n_ism*f_nh; // cm^-3 -- 30Dor
//const double nhe = n_ism*f_nhe; // cm^-3 -- 30Dor

double n_250 = n_ism*1./(250);
double nh_neut = nh*f_nh_neut; // cm^-3
double nh_ion = nh*f_nh_ion; // cm^-3
double nhi = nh_neut*f_nh_neut_atom; // cm^-3
double nh_mol = nh_neut*f_nh_neut_mol*1./2.; // cm^-3
//const double nhi = gas_density*f_nh*f_nh_neut*f_atom_h*1./mass_hyd_grams; // cm^-3
//const double nh_mol = gas_density*f_nh*f_nh_neut*f_mol_h*1./(2.*mass_hyd_grams); // cm^-3

//const double snrate = 0.08/year; // s^-1 -- NGC 253
double snrate = 0.02/year; // s^-1 -- NGC 253 -- Lacki et al. 2012
//const double snrate = 0.15/year; // s^-1 -- 30Dor

double snenergy = pow(10.0,51.0)*ergs2eV; // eV
//const double sneff_proton = 0.1;
//const double sneff_proton = 0.2; // Lacki et al. 2012 leptonic
double sneff_proton = 0.1; // Lacki et al. 2012 hadronic
//const double nh = 1;

/*// Proton injection spectrum parameters - power law
const double prot_p = 2.1;//2.1; 2.75;  // cf. Chakraborty & Fields 2013 for starburst galaxies -- 2.1; C&F say that range is: 2.0--2.4; Lacki et al. say the range is 2.0--2.6
//const double elec_norm_pl = ;
// 30Dor range: 2.1--2.8 (Murphy et al.)

// Proton injection spectrum parameters - smoothly broken power law
const double prot_p_faint = 1.8;     // cf. Strong et al. 2010
const double prot_p_bright = 2.25;
const double prot_n = 1.0; // strength of the transition
const double prot_break = 4.e9; // eV
const double prot_high_cutoff = 1.e15; // eV
//const double elec_norm_sbpl = ;*/

// Range of Proton injection spectrum
const double prot_E_min = protonmass; // eV
const double prot_E_max = /*(1.e6)*protonmass*/1.e20; // eV
const double log_prot_gmin = log10(prot_E_min*1./protonmass);
const double log_prot_gmax = log10(prot_E_max*1./protonmass);
const double zcr = 1.0;

// Starting proton spectrum -- dN/dVd(log(gamma))

TGraph *prdist = new TGraph("/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_exp_cutoff_w_args.dat");
//TGraph *prdist = new TGraph("Tests/Lacki/50muG/ne_prot_ss_finite_exp_cutoff_w_args.dat");

const double pi0_E_max = K_pi*(prot_E_max - protonmass); // eV
const double log_E_pi0_max = log10(pi0_E_max);
const double gam_pi0_max = pi0_E_max*1./pi0mass;
const double log_gam_pi0_max = log10(gam_pi0_max);
const double pion_E_max = pi0_E_max;
const double log_E_pion_max = log10(pion_E_max);
const double gam_pion_max = pion_E_max*1./pionmass;
const double log_gam_pion_max = log10(gam_pion_max);
const double gam_muon_max = pion_E_max*1./muonmass;
const double log_gam_muon_max = log10(gam_muon_max);

// Range of electron injection spectrum
//const double elec_E_min = 1.5e6; // eV
const double elec_E_min = elecmass;
const double elec_E_max = /*(1.e6)*elecmass*/1.e17; // eV
//const double elec_E_max = 1.e17; // eV
const double log_elec_gmin = log10(elec_E_min*1./elecmass);
const double log_elec_gmax = log10(elec_E_max*1./elecmass);

//const double specflag = 2.0; // Must be 1 or 2 at the moment

// Photon Backgrounds
// Note that all differentials that have to be integrated over energy are to be integrated over log_10(energy) rather than energy. So, if they look a little strange, that would be the reason.

const double isrfmodeldistmpc = 50; // NGC 253, but you don't need to worry about it here

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
double f_e_high(double x_e, double protenergy) // x_e = E_e/E_p; [protenergy] = eV; f_nu_e_high ~ f_e_high
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
double f_nu_mu_high(double x_nu_mu, double protenergy) // x_nu_mu = E_nu_mu/E_p; [protenergy] = eV
{
    double prot_energy_tev = protenergy*(1.e-12); // TeV
    double ln_p_energy_tev = log(prot_energy_tev);
    double ln_x_nu_mu = log(x_nu_mu);
    double y_nu_mu = x_nu_mu*1./0.427;
    double ln_y_nu_mu = log(y_nu_mu);
    double B_prime, beta_prime, ynu_pow_betapr, k_prime;
    double factor1, factor1_pow_4, factor2, factor3;
    double F_nu_1, F_nu_2;
    
    B_prime = 1.75 + 0.204*ln_p_energy_tev + 0.010*ln_p_energy_tev*ln_p_energy_tev;
    beta_prime = 1./(1.67 + 0.111*ln_p_energy_tev+0.0038*ln_p_energy_tev*ln_p_energy_tev);
    k_prime = 1.07 - 0.086*ln_p_energy_tev + 0.002*ln_p_energy_tev*ln_p_energy_tev;
    
    ynu_pow_betapr = pow(y_nu_mu,beta_prime);
    
    factor1 = (1.-ynu_pow_betapr)*1./(1.+k_prime*ynu_pow_betapr*(1.-ynu_pow_betapr));
    factor1_pow_4 = pow(factor1,4.0);
    factor2 = 1./ln_y_nu_mu - 4.*beta_prime*ynu_pow_betapr/(1.-ynu_pow_betapr);
    factor3 = 4.*k_prime*beta_prime*ynu_pow_betapr*(1.-2.*ynu_pow_betapr)/(1.+k_prime*ynu_pow_betapr*(1.-ynu_pow_betapr));
    
    if (x_nu_mu < 0.427) F_nu_1 = B_prime*(ln_y_nu_mu*1./y_nu_mu)*factor1_pow_4*(factor2 - factor3);
    else F_nu_1 = 0;
    //F_nu_1 = B_prime*(ln_y_nu_mu*1./y_nu_mu)*factor1_pow_4*(factor2 - factor3);
    F_nu_2 = f_e_high(x_nu_mu,protenergy);
    
    return F_nu_1 + F_nu_2;
}
double f_nu_e_high(double x_nu_e, double protenergy) // x_nu_e = E_nu_e/E_p; [protenergy] = eV
{
    return f_e_high(x_nu_e,protenergy);
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
    
    return factor1*factor2fourth*factor3*factor4;
}
double f_eta_high(double x_eta, double protenergy) // x_eta = E_eta/E_p; [protenergy] = eV
{
    double ln_x = log(x_eta);
    
    return (0.55 + 0.028*ln_x)*(1.-etamass*1./(x_eta*protenergy))*f_pi_high(x_eta,protenergy);
}
double e_f_pi_high_diff(double *x_pi, double *par)
{
    double protenergy = par[0];
    double E_pi = x_pi[0]*protenergy;
    
    return E_pi*f_pi_high(x_pi[0],protenergy);
}
double e_f_eta_high_diff(double *x_eta, double *par)
{
    double protenergy = par[0];
    double E_eta = x_eta[0]*protenergy;
    
    return E_eta*f_eta_high(x_eta[0],protenergy);
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

double kn_elec_spectrum_tor(double gamma_prim, double gamma_e) // Currently not in use; for testing purposes only
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
 f_of_t_hi = 1. - prim_beta_sqrd*(T_e*1./T_max) + 0.5*(T_e*1./prim_energy)*(T_e*1./prim_energy); // spin 1/2 formula -- other formulas are given in the notes and in High Energy Particles by Bruno Rossi (1952) or the Particle Data Group Review
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

// Proton Distribution

double NE_prot(double loggamma)
{
    //TGraph *prdist = new TGraph("Diffusion_Only_Models/10Myr/ne_prot_ss_finite_with_he.dat"); // dN/dVd(log(gamma))
    
    //double protenergy_tev = pow(10.0,loggamma)*protonmass*(1.e-12);
    
    // -- Using these for testing --
    
    //return Jp_kel(protenergy_tev)*(protenergy_tev*log(10.0)); // dN/dVd(log(gamma))
    
    // -- Using either of the these two for results --
    
    return prdist->Eval(loggamma); // Might try spline
    //return prdist->Eval(loggamma,0,"S");
    
    // -- Not using the stuff below --
    
    /*ifstream fpdist("Diffusion_Only_Models/10Myr/ne_prot_ss_finite_with_he.dat"); // dN/dVd(log(gamma))
    
    vector<double> logdistgamma, dNdloggamdv;
    int values = 0, index = 0, i;
    double gamvalue, distvalue, p, b;
    bool doneflag = false;
    
    //double logdistgamma[lines], dNdloggamdv[lines];
    while (!fpdist.eof()) {
        fpdist >> gamvalue >> distvalue;
        //cout << gamvalue << " " << distvalue << endl;
        logdistgamma.push_back(gamvalue);
        dNdloggamdv.push_back(distvalue);
        values++;
    }
    
    values--;
    
    fpdist.close();
    
    //if (loggamma < 0.1) cout << loggamma << endl;
    
    if ((loggamma >= logdistgamma[0])&&(loggamma <= logdistgamma[values-1])){
        for (i = 0; i < values; i++) {
            if ((loggamma >= logdistgamma[i])&&(loggamma < logdistgamma[i+1])) {
                index = i;
                if (loggamma == logdistgamma[i]) doneflag = true;
                break;
            }
        }
        if (doneflag) return dNdloggamdv[index];
        else {
            p = (loggamma - logdistgamma[index])*1./(logdistgamma[index+1]-logdistgamma[index]);
            return (1.-p)*dNdloggamdv[index]+p*dNdloggamdv[index+1];
        }
    }
    else{
        //cout << loggamma << endl;
        if (loggamma < logdistgamma[0]) {
            p = (dNdloggamdv[1]-dNdloggamdv[0])*1./(logdistgamma[1]-logdistgamma[0]);
            b = dNdloggamdv[0] - p*logdistgamma[0];
            
            if (loggamma < 0.1) cout << "Here: " << loggamma << " " << p*loggamma + b << endl;
            
            return p*loggamma + b;
        }
        if (loggamma > logdistgamma[values-1]) {
            p = (dNdloggamdv[values-1]-dNdloggamdv[index-2])*1./(logdistgamma[values-1]-logdistgamma[values-2]);
            b = dNdloggamdv[values-1] - p*logdistgamma[values-1];
            //cout << loggamma << " " << p*loggamma + b << endl;
            return p*loggamma + b;
        }
    }*/
}

// Total Secondary Particle Spectra Differentials

double q_pi0_high_diff(double *x_gamma, double *par) // x_gamma = E_gamma/E_p;
{
    double gammaenergy = par[0]; // eV
    double protenergy = gammaenergy*1./x_gamma[0];
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    
    if (prot_gamma > 1) return pp_inel_high(protenergy)*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*f_gamma_high(x_gamma[0],protenergy)*(1./x_gamma[0]); // photons cm^-1 eV^-1
    //if (prot_gamma > 1) return pp_inel_high(protenergy)*Jp*f_gamma_high(x_gamma[0],protenergy)*(1./x_gamma[0]); // photons cm^-1 eV^-1
    else return 0;
}
double q_pion_high_diff(double *x_e, double *par) // x_e = E_e/E_p
{
    double elecenergy = par[0]; // eV
    double protenergy = elecenergy*1./x_e[0];
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    
    //cout << x_e[0] << " " << logprotgamma << " " << NE_prot(logprotgamma) << endl;
    
    if (prot_gamma > 1) return pp_inel_high(protenergy)*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*f_e_high(x_e[0],protenergy)*(1./x_e[0]); // electrons cm^-1 eV^-1
    //if (prot_gamma > 1) return pp_inel_high(protenergy)*Jp*f_e_high(x_e[0],protenergy)*(1./x_e[0]); // electrons cm^-1 eV^-1
    else return 0;
    
}
double q_gamma_high(double gamenergy) // eV
{
    double qpiinteg;
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 q_pi_h_diff("q_pi_h_diff",q_pi0_high_diff,-0.5,1.5,1);
    q_pi_h_diff.SetParameter(0,gamenergy);
    q_pi_h_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_pi_h_diff.IntegralFast(np,x,w,(1.0e-3),1.0); // photons cm^-1 eV^-1
    
    //cout << qpiinteg << endl;
    
    delete [] x;
	delete [] w;
    
    return qpiinteg*speedoflight*nh; // photons cm^-3 s^-1 eV^-1
}
double q_e_high(double elecenergy) // eV
{
    double qpiinteg;
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 q_pi_h_diff("q_pi_h_diff",q_pion_high_diff,-0.5,1.5,1);
    q_pi_h_diff.SetParameter(0,elecenergy);
    q_pi_h_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_pi_h_diff.IntegralFast(np,x,w,(1.0e-3),1.0); // electrons cm^-1 eV^-1
    
    delete [] x;
	delete [] w;
    
    return qpiinteg*speedoflight*nh; // electrons cm^-3 s^-1 eV^-1
}
double q_nu_mu_high_diff(double *x_nu_mu, double *par) // x_nu_mu = E_nu_mu/E_p
{
    double numuenergy = par[0]; // eV
    double protenergy = numuenergy*1./x_nu_mu[0];
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    
    if (prot_gamma > 1) return pp_inel_high(protenergy)*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*f_nu_mu_high(x_nu_mu[0],protenergy)*(1./x_nu_mu[0]); // neutrinos cm^-1 eV^-1
    //if (prot_gamma > 1) return pp_inel_high(protenergy)*Jp*f_nu_mu_high(x_nu_mu[0],protenergy)*(1./x_nu_mu[0]); // neutrinos cm^-1 eV^-1
    else return 0;
}
double q_nu_mu_high(double numuenergy) // eV
{
    double qpiinteg;
    
    int np = 50;
    double *x = new double[np];
    double *w = new double[np];
    
    TF1 q_nm_h_diff("q_nm_h_diff",q_nu_mu_high_diff,-0.5,1.5,1);
    q_nm_h_diff.SetParameter(0,numuenergy);
    q_nm_h_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_nm_h_diff.IntegralFast(np,x,w,(1.0e-3),1.0);
    
    delete [] x;
	delete [] w;
    
    return qpiinteg*speedoflight*nh; // neutrinos cm^-3 s^-1 eV^-1
}
double q_nu_e_high_diff(double *x_nu_e, double *par) // x_nu_e = E_nu_e/E_p
{
    double nueenergy = par[0]; // eV
    double protenergy = nueenergy*1./x_nu_e[0];
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    
    if (prot_gamma > 1) return pp_inel_high(protenergy)*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*f_nu_e_high(x_nu_e[0],protenergy)*(1./x_nu_e[0]); // neutrinos cm^-1 eV^-1
    //if (prot_gamma > 1) return pp_inel_high(protenergy)*Jp*f_nu_e_high(x_nu_e[0],protenergy)*(1./x_nu_e[0]); // neutrinos cm^-1 eV^-1
    else return 0;
}
double q_nu_e_high(double nueenergy) // eV
{
    double qpiinteg;
    
    int np = 50;
    double *x = new double[np];
    double *w = new double[np];
    
    TF1 q_ne_h_diff("q_ne_h_diff",q_nu_e_high_diff,-0.5,1.5,1);
    q_ne_h_diff.SetParameter(0,nueenergy);
    q_ne_h_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_ne_h_diff.IntegralFast(np,x,w,(1.0e-3),1.0);
    
    delete [] x;
	delete [] w;
    
    return qpiinteg*speedoflight*nh; // neutrinos cm^-3 s^-1 eV^-1
}
double q_pi0_low_diff(double *logEpi0, double *par)
{
    double E_pi0 = pow(10.0,logEpi0[0]);
    double protenergy = protonmass + E_pi0*1./K_pi;
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    double denom = sqrt(E_pi0*E_pi0 - pi0mass*pi0mass);
    
    //return q_gamma_p_pi0_low(E_pi0)*Jp*(1./denom)*E_pi0*log(10.0); // pions cm^-3 s^-1 eV^-1 (log(eV))^-1
    
    return q_gamma_p_pi0_low(E_pi0)*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*(1./denom)*E_pi0*log(10.0); // pions cm^-3 s^-1 eV^-1 (log(eV))^-1
}
double q_pion_low_diff(double *loggampion, double *par)
{
    double gampion = pow(10.0,loggampion[0]);
    double E_pion = gampion*pionmass; // eV
    double protenergy = protonmass + E_pion*1./K_pi;
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    double denom = sqrt(gampion*gampion - 1.);
    
    //return q_gamma_p_pi0_low(E_pion)*pionmass*Jp*(1./denom)*gampion*log(10.0); // pions cm^-3 s^-1 (log(gamma))^-1
    
    return q_gamma_p_pi0_low(E_pion)*pionmass*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*(1./denom)*gampion*log(10.0); // pions cm^-3 s^-1 (log(gamma))^-1
}
double q_muon_low_diff(double *loggammuon, double *par)
{
    double gammuon = pow(10.0,loggammuon[0]);
    double E_pion = gammuon*pionmass; // eV
    double protenergy = protonmass + E_pion*1./K_pi;
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    double denom = sqrt(gammuon*gammuon - 1.);
    
    //return q_gamma_p_pi0_low(E_pion)*pionmass*Jp*(1./denom)*gammuon*log(10.0); // muons cm^-3 s^-1 gamma^-1 (log(gamma))^-1
    
    return q_gamma_p_pi0_low(E_pion)*pionmass*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*(1./denom)*gammuon*log(10.0); // muons cm^-3 s^-1 gamma^-1 (log(gamma))^-1
}
double q_pion_no_denom_low_diff(double gampion)
{
    //double gampion = pow(10.0,loggampion);
    double E_pion = gampion*pionmass; // eV
    double protenergy = protonmass + E_pion*1./K_pi;
    double prot_gamma = protenergy*1./protonmass;
    double logprotgamma = log10(prot_gamma);
    
    //return q_gamma_p_pi0_low(E_pion)*pionmass*Jp*gampion*log(10.0); // pions cm^-3 s^-1 (log(gamma))^-1
    
    return q_gamma_p_pi0_low(E_pion)*pionmass*NE_prot(logprotgamma)*(1./(protenergy*log(10.0)))*gampion*log(10.0); // pions cm^-3 s^-1 (log(gamma))^-1
}
double q_gamma_low(double gamenergy) // eV
{
    double E_min = gamenergy + pi0mass*pi0mass*1./(4.*gamenergy);
    double logEmin = log10(E_min);
    double qpiinteg;
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    //cout << E_min << endl;
    
    TF1 q_pi_l_diff("q_pi_l_diff",q_pi0_low_diff,0.0,20.0,1);
    q_pi_l_diff.SetParameter(0,0);
    q_pi_l_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_pi_l_diff.IntegralFast(np,x,w,logEmin,log_E_pi0_max); // pions cm^-3 s^-1 eV^-1
    
    //cout << logEmin << " " << log_pi0_gmax << endl;
    
    delete [] x;
	delete [] w;
    
    return 2.0*qpiinteg; // photons cm^-3 s^-1 eV^-1 -- Note: Still must be normalized
}
double q_pi_nu_mu_low(double numuenergy) // eV
{
    double numuenergy_cms = (pionmass*pionmass - muonmass*muonmass)*1./(2.0*pionmass);
    double gampi_min = (1./2.)*((numuenergy_cms*1./numuenergy) + (numuenergy*1./numuenergy_cms));
    double loggampimin = log10(gampi_min);
    double qpiinteg;
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 q_pi_l_diff("q_pi_l_diff",q_pion_low_diff,0.0,20.0,1);
    q_pi_l_diff.SetParameter(0,0);
    q_pi_l_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_pi_l_diff.IntegralFast(np,x,w,loggampimin,log_gam_pion_max); // pions cm^-3 s^-1
    
    delete [] x;
	delete [] w;
    
    return 2.0*(1./(2.*numuenergy_cms))*qpiinteg; // neutrinos cm^-3 s^-1 eV^-1 -- Note: (1) Still must be normalized; (2) includes both muon neutrinos and muon antineutrinos from the pion decays
}
double q_muon_chi_low_diff(double *chi, double *par)
{
    //double loggammuon = par[0];
    //double gammuon = pow(10.0,loggammuon);
    double gammuon = par[0];
    double eps1 = 3.68, eps2 = 3.55*gammuon;
    double gampion = (1./(eps1*eps1-chi[0]*chi[0]))*(eps2*eps1+chi[0]*sqrt(eps2*eps2-eps1*eps1+chi[0]*chi[0]));
    double dgampi_dgammu = (1./(eps1*eps1-chi[0]*chi[0]))*(3.55*eps1 + chi[0]*(1./sqrt(eps2*eps2-eps1*eps1+chi[0]*chi[0]))*(3.55*eps2));
    double dloggampi_dloggammu = (gammuon*1./gampion)*dgampi_dgammu;
    
    return q_pion_no_denom_low_diff(gampion)*dloggampi_dloggammu; // pions cm^-3 s^-1 (log(gamma))^-1
}
double f_dist_nu(double nuenergy, double gammuon, int nuflag)
{
    double zeta = nuenergy*1./(gammuon*muonmass);
    double zeta_sqrd = zeta*zeta;
    double beta_mu = sqrt(1.-1./(gammuon*gammuon));
    double beta_mu_sqrd = beta_mu*beta_mu;
    double beta_mu_pl1_pow_3 = pow((1.+beta_mu),3.0);
    double gammuon_sqrd = gammuon*gammuon;
    double gammuon_pow_5 = pow(gammuon,5.0);
    double f_nu = 0;
    
    // nuflag = 1 for muon and electron neutrinos
    // nuflag = 2 for muon and electron antineutrinos
    
    switch (nuflag) {
        case 1:
            if ((zeta >= 0)&&(zeta <= ((1.-beta_mu)/2.))){
                f_nu = 16.*gammuon_pow_5*zeta_sqrd*(1./muonmass)*(3./gammuon_sqrd-(4./3.)*zeta*(3.+beta_mu_sqrd));
            }
            if ((zeta > ((1.-beta_mu)/2.))&&(zeta <= ((1.+beta_mu)/2.))){
                f_nu = (1./(beta_mu*gammuon*muonmass))*((5./3.)+(4./beta_mu_pl1_pow_3)*zeta_sqrd*(8.*zeta/3. - 3.*(1.+beta_mu)));
            }
            break;
        case 2:
            if ((zeta >= 0)&&(zeta <= ((1.-beta_mu)/2.))){
                f_nu = 32.*gammuon_pow_5*zeta_sqrd*(1./muonmass)*(3./gammuon_sqrd-2.*zeta*(3.+beta_mu_sqrd));
            }
            if ((zeta > ((1.-beta_mu)/2.))&&(zeta <= ((1.+beta_mu)/2.))){
                f_nu = (1./(2.*beta_mu*gammuon*muonmass))*(1.+(4./beta_mu_pl1_pow_3)*zeta_sqrd*(4.*zeta - 3.*(1.+beta_mu)));
            }
    }
    return f_nu;
}
double q_muon_integ_chi_low_diff(double *loggammuon, double *par)
{
    double nuenergy = par[0];
    int nuflag = (int)(par[1]+0.5);
    
    double gammuon = pow(10.0,loggammuon[0]);
    double f_nu = f_dist_nu(nuenergy,gammuon,nuflag);
    double eps1 = 3.68, eps2 = 3.55*gammuon;
    double chi_min, qpiinteg;
    
    if (eps1 <= eps2) chi_min = -1.0;
    else chi_min = sqrt(eps1*eps1-eps2*eps2);
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 q_mu_chi_l_diff("q_mu_chi_l_diff",q_muon_chi_low_diff,-1.5,1.5,1);
    q_mu_chi_l_diff.SetParameter(0,gammuon);
    q_mu_chi_l_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qpiinteg = q_mu_chi_l_diff.IntegralFast(np,x,w,chi_min,1.0);
    
    delete [] x;
	delete [] w;
    
    return f_nu*qpiinteg; // neutrinos cm^-3 s^-1 eV^-1 (log(gamma))^-1
}
double q_mu_dec_nu_low(double nuenergy, double nuflag_db) // eV; Note that nuflag_db must be equal to 1 or 2
{
    double gammuon_min = (nuenergy*1./muonmass)+(muonmass*1./(4.*nuenergy));
    double log_gammu_min = log10(gammuon_min);
    double qinteg;
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 q_mu_integ_chi_l_diff("q_mu_integ_chi_l_diff",q_muon_integ_chi_low_diff,0.0,20.0,2);
    q_mu_integ_chi_l_diff.SetParameters(nuenergy,nuflag_db);
    q_mu_integ_chi_l_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qinteg = q_mu_integ_chi_l_diff.IntegralFast(np,x,w,log_gammu_min,log_gam_muon_max);
    
    delete [] x;
	delete [] w;
    
    return 0.5*qinteg; // neutrinos cm^-3 s^-1 eV^-1
}
double q_nu_mu_low(double numuenergy)
{
    return (q_pi_nu_mu_low(numuenergy) + q_mu_dec_nu_low(numuenergy,1.0) + q_mu_dec_nu_low(numuenergy,2.0)); // neutrinos cm^-3 s^-1 eV^-1-- Note: (1) Still must be normalized; (2) includes both muon neutrinos and muon antineutrinos from the pion decays
}
double q_nu_e_low(double nueenergy)
{
    return (q_mu_dec_nu_low(nueenergy,1.0) + q_mu_dec_nu_low(nueenergy,2.0)); // neutrinos cm^-3 s^-1 eV^-1 -- Note: (1) Still must be normalized; (2) includes both muon neutrinos and muon antineutrinos from the pion decays
}
double q_e_low_diff(double *loggamepr, double *par)
{
    double gamma_e_prime = pow(10.0,loggamepr[0]);
    double gam_e_pr_sqrd = gamma_e_prime*gamma_e_prime;
    double gamma_e = par[0];
    double gamma_e_sqrd = gamma_e*gamma_e;
    double gamma_mu_1, gamma_mu_2, gamma_e_pr_max = 104., p_game_prime;
    double gamma_e_pr_max_cube = gamma_e_pr_max*gamma_e_pr_max*gamma_e_pr_max;
    double loggam_mu_1, loggam_mu_2;
    double denom = sqrt(gam_e_pr_sqrd - 1.);
    double qmuinteg;
    
    gamma_mu_1 = gamma_e*gamma_e_prime - sqrt(gamma_e_sqrd - 1.)*sqrt(gam_e_pr_sqrd - 1.);
    gamma_mu_2 = gamma_e*gamma_e_prime + sqrt(gamma_e_sqrd - 1.)*sqrt(gam_e_pr_sqrd - 1.);
    
    loggam_mu_1 = log10(gamma_mu_1);
    loggam_mu_2 = log10(gamma_mu_2);
    
    p_game_prime = 2.*gam_e_pr_sqrd*(3.-2.*gamma_e_prime/gamma_e_pr_max)*(1./(gamma_e_pr_max_cube));
    
    int np = 50;
	double *x = new double[np];
	double *w = new double[np];
    
    TF1 q_mu_diff("q_mu_diff",q_muon_low_diff,loggam_mu_1-0.5,loggam_mu_2+0.5,1);
    q_mu_diff.SetParameter(0,0.0);
    q_mu_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qmuinteg = q_mu_diff.IntegralFast(np,x,w,loggam_mu_1,loggam_mu_2); // pions cm^-3 s^-1 gamma^-1
    
    delete [] x;
	delete [] w;
    
    //if (qmuinteg <= 0) cout << loggam_mu_1 << " " << loggam_mu_2 << endl;
    
    //cout << qmuinteg << endl;
    
    return 0.5*p_game_prime*qmuinteg*gamma_e_prime*log(10.0)/denom; // pions cm^-3 s^-1 gamma^-1 (log(gamma))^-1
}
double q_e_low(double elecenergy) // eV
{
    double gamma_e = elecenergy*1./elecmass;
    double gamma_e_pr_max = 104;
    double loggam_e_pr_max = log10(gamma_e_pr_max);
    double qeinteg;
    
    int np = 50;
    double *x = new double[np];
    double *w = new double[np];
    
    //cout << E_min << endl;
    
    TF1 q_e_l_diff("q_e_l_diff",q_e_low_diff,-0.5,3.0,1);
    q_e_l_diff.SetParameter(0,gamma_e);
    q_e_l_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qeinteg = q_e_l_diff.IntegralFast(np,x,w,0.0,loggam_e_pr_max); // pions cm^-3 s^-1 gamma^-1
    
    delete [] x;
    delete [] w;
    
    return 2.0*qeinteg/elecmass; // electrons cm^-3 s^-1 eV^-1 -- Note: Still must be normalized
}
double q_gamma(double gamenergy) // eV
{
    double normalization;
    double qhigh, qlow;
    
    normalization = q_gamma_high(1.e11)*1./(q_gamma_low(1.e11));
    
    //cout << normalization << endl;
    
    qhigh = q_gamma_high(gamenergy);
    qlow = normalization*q_gamma_low(gamenergy);
    
    //cout << normalization << endl;
    
    //cout << gamenergy << " " << qlow << " " << qhigh << endl;
    
    if (gamenergy <= (1.e11)) return qlow;
    else return qhigh;
    
    //return (qhigh + qlow); // photons cm^-3 s^-1 eV^-1
}
double q_nu_mu(double numuenergy) // eV
{
    double normalization;
    
    normalization = q_nu_mu_high(1.e11)*1./(q_nu_mu_low(1.e11));
    if (numuenergy <= 1.e11) return normalization*q_nu_mu_low(numuenergy);
    else return q_nu_mu_high(numuenergy); // neutrinos cm^-3 s^-1 eV^-1
}
double q_nu_e(double nueenergy) // eV
{
    double normalization;
    
    normalization = q_nu_e_high(1.e11)*1./(q_nu_e_low(1.e11));
    if (nueenergy <= 1.e11) return normalization*q_nu_e_low(nueenergy);
    else return q_nu_e_high(nueenergy); // neutrinos cm^-3 s^-1 eV^-1
}
double q_e_pp(double elecenergy) // eV
{
    double normalization;
    
    normalization = q_e_high(1.e11)*1./(q_e_low(1.e11));
    
    //cout << normalization << endl;
    
    //cout << q_e_low(elecenergy) << " " << q_e_high(elecenergy) << endl;
    
    if (elecenergy <= 1.e11) return normalization*q_e_low(elecenergy);
    else return q_e_high(elecenergy); // electrons s^-1 cm^-3 eV^-1
}
double q_e_kn_diff(double *loggamma_p, double *par)
{
    double gamma_p = pow(10.0,loggamma_p[0]);
    double protenergy = gamma_p*protonmass;
    double gamma_e = par[0];
    
    return kn_elec_spectrum_pdg(gamma_p,gamma_e)*NE_prot(loggamma_p[0])*(1./(protenergy*log(10.0)))*protenergy*log(10.0)*1./elecmass; // dN/dVdtdE_ed(log(gamma_p)) -- electrons s^-1 eV^-1 cm^-3
}
double q_e_kn(double elecenergy) // eV
{
    double gamme_e = elecenergy*1./elecmass;
    double qkninteg;
    double gamma_p1, gamma_p2, gamma_p_min, loggamma_p_min;
    double s = elecmass*1./(atom_mass_hyd*protonmass); // if using higher Z primaries, must replace protonmass with the mass of the primary
    
    gamma_p1 = 0.5*((gamme_e-1.)*s+sqrt((gamme_e-1.)*(gamme_e-1.)*s*s+4*(0.5*(gamme_e-1.)*(1.+s*s)-1.)));
    gamma_p2 = 0.5*((gamme_e-1.)*s-sqrt((gamme_e-1.)*(gamme_e-1.)*s*s+4*(0.5*(gamme_e-1.)*(1.+s*s)-1.)));
    
    if (((gamme_e-1.)*(gamme_e-1.)*s*s + 4*0.5*(gamme_e-1.)*(1.+s*s)) > 4) {
        if (gamma_p1 > 0.0) gamma_p_min = gamma_p1;
        else gamma_p_min = gamma_p2;
        if (log10(gamma_p_min) > 0) loggamma_p_min = log10(gamma_p_min);
        else loggamma_p_min = log_prot_gmin;
        //cout << gamma_p1 << " " << gamma_p2 << endl;
    }
    else loggamma_p_min = log_prot_gmin;
    
    //loggamma_p_min = log10(gamma_p_min);
    
    //cout << loggamma_p_min << endl;
    
    int np = 50;
    double *x = new double[np];
    double *w = new double[np];
    
    TF1 q_kn_diff("q_kn_diff",q_e_kn_diff,(log_prot_gmin-0.5),(log_prot_gmax+0.5),1);
    q_kn_diff.SetParameter(0,gamme_e);
    q_kn_diff.CalcGaussLegendreSamplingPoints(np,x,w,1.e-15);
    qkninteg = q_kn_diff.IntegralFast(np,x,w,loggamma_p_min,log_prot_gmax); // dN/dVdtdE_e -- electrons s^-1 cm^-3 eV^-1
    
    //cout << loggamma_p_min << " " << qkninteg << endl;
    
    delete [] x;
    delete [] w;
    
    return qkninteg;
}
double q_elec(double elecenergy) // eV
{
    //cout << q_e_pp(elecenergy) << " " << q_e_kn(elecenergy) << endl;
    
    return (0.5*q_e_pp(elecenergy)+q_e_kn(elecenergy)); // dN/dVdtdE_e -- electrons s^-1 cm^-3 eV^-1
    //return 0.5*q_e_pp(elecenergy);
}
double q_posi(double elecenergy) // eV
{
    return 0.5*q_e_pp(elecenergy); // dN/dVdtdE_e -- electrons s^-1 cm^-3 eV^-1
}
int main(){
//int main(int argc, char **argv){
    
    // Parameters to be input on the command line
    // Emission flag: [1.0] Hadronic; [2.0] Leptonic
    
    // To run the current model, type: |> ./<name of executable>.out 1.0
    
    /*ofstream fne("Lacki/400muG/qe_e_minus_ss_finite_exp_cutoff.dat");
    ofstream fposi("Lacki/400muG/qe_e_plus_ss_finite_exp_cutoff.dat");
    ofstream fgam("Lacki/400muG/pp_gam_spec_ss_finite_exp_cutoff.dat");*/
    
     
    
    
    double emispar = 1.0; //atof(argv[1]);
    
    int emisflag = (int)(emispar+0.5);
    
    char fnename[250];
    char fposiname[250];
    char fgamname[250];
    char fnuename[250];
    char fnumuname[250];
    
    switch (emisflag) {
        case 1: // NGC 253 -- Lacki et al. 2012 hadronic
            irscaleheight = 5.0*pc; //[cm]
            regBfield = 4.e-4; //[G]
            //ofstream fne("Tests/Lacki/400muG/qe_e_minus_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fnename,"Tests/Lacki/400muG/qe_e_minus_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fposi("Tests/Lacki/400muG/qe_e_plus_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fposiname,"Tests/Lacki/400muG/qe_e_plus_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fgam("Tests/Lacki/400muG/pp_gam_spec_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fgamname,"Tests/Lacki/400muG/pp_gam_spec_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fnue("Tests/Lacki/400muG/qe_nu_e_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fnuename,"Tests/Lacki/400muG/qe_nu_e_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fnumu("Tests/Lacki/400muG/qe_nu_mu_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fnumuname,"Tests/Lacki/400muG/qe_nu_mu_ss_finite_exp_cutoff_w_args.dat");
            //prdist = new TGraph("Tests/Lacki/400muG/ne_prot_ss_finite_exp_cutoff_w_args.dat");
            break;
        case 2: // NGC 253 -- Lacki et al. 2012 leptonic
            irscaleheight = 50.0*pc;
            regBfield = 5.e-5; //[G]
            //ofstream fne("Tests/Lacki/50muG/qe_e_minus_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fnename,"Tests/Lacki/50muG/qe_e_minus_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fposi("Tests/Lacki/50muG/qe_e_plus_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fposiname,"Tests/Lacki/50muG/qe_e_plus_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fgam("Tests/Lacki/50muG/pp_gam_spec_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fgamname,"Tests/Lacki/50muG/pp_gam_spec_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fnue("Tests/Lacki/50muG/qe_nu_e_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fnuename,"Tests/Lacki/50muG/qe_nu_e_ss_finite_exp_cutoff_w_args.dat");
            //ofstream fnumu("Tests/Lacki/50muG/qe_nu_mu_ss_finite_exp_cutoff_w_args.dat");
            //sprintf(fnumuname,"Tests/Lacki/50muG/qe_nu_mu_ss_finite_exp_cutoff_w_args.dat");
            //prdist = new TGraph("Tests/Lacki/50muG/ne_prot_ss_finite_exp_cutoff_w_args.dat");
            break;
    }
    
    /*ofstream fne(fnename);
    ofstream fposi(fposiname);
    ofstream fgam(fgamname);
    ofstream fnue(fnuename);
    ofstream fnumu(fnumuname);*/
    
    ofstream fne;
    ofstream fposi;
    ofstream fgam;
    ofstream fnue;
    ofstream fnumu; 
    
    //double energy, gamma, loggamma;
    double logbinsize = .1;//, loggammamin = 1.0, loggammamax = 12.0;
    //double loggamma = loggammamin;
    double loggamma_e = log_elec_gmin+logbinsize;
    double qe_e_ss, qe_p_ss, qe_nue_ss, qe_numu_ss, normalization, qekininteg, qe_gamma;
    double logepsmin = 6.0, logepsbinsize = .1, logepsmax = 15.0;
    double logeps = logepsmin, eps, E_elec;
    double inputplleindex = 1.9, inputlogEcut = pow(10.0,15.0);
    //double modelpars[7];
    
    char protdistfilename[250];
    
    //double emispar = atof(argv[1]), specpar = atof(argv[2]), inputplleindex = atof(argv[3]), inputbplhiindex = atof(argv[4]), inputbpln = atof(argv[5]), inputbplEbreak = pow(10.0,atof(argv[6])), inputlogEcut = pow(10.0,atof(argv[7]));
    
    // Resetting certain source paramaters
    
    h_100 = irscaleheight*1./(100.*pc);
    volirregion = 2.*pi*irregion*irregion*irscaleheight; // NGC 253 -- Disk volume
    gas_density = gas_mass*1./volirregion; // [g cm^-3]
    nh = gas_density*f_nh*1./mass_hyd_grams; // [cm^-3]
    nhe = gas_density*f_nhe*1./mass_he_grams; // [cm^-3]
    n_ism = nh + nhe; // [cm^-3]
    n_250 = n_ism*1./(250);
    nh_neut = nh*f_nh_neut; // [cm^-3]
    nh_ion = nh*f_nh_ion; // [cm^-3]
    nhi = nh_neut*f_nh_neut_atom; // [cm^-3]
    nh_mol = nh_neut*f_nh_neut_mol*1./2.; // [cm^-3]
    
    
    
    // inputplleindex loop starts here
    
    for (inputplleindex=2.55; inputplleindex<=2.625; inputplleindex+=0.025)
    {  
    
    
    
        sprintf(protdistfilename,"/Users/CEE/Desktop/NGC_253/ne_prot_ss_finite_B_%.1f_plindex_%.3f_Ecut_%.1f.dat",(regBfield*1.e6),inputplleindex,log10(inputlogEcut));
        prdist = new TGraph(protdistfilename);
    
        sprintf(fnename,"/Users/CEE/Desktop/NGC_253/qe_e_minus_ss_finite_B_%.1f_plindex_%.3f_Ecut_%.1f.dat",(regBfield*1.e6),inputplleindex,log10(inputlogEcut));
        sprintf(fposiname,"/Users/CEE/Desktop/NGC_253/qe_e_plus_ss_finite_B_%.1f_plindex_%.3f_Ecut_%.1f.dat",(regBfield*1.e6),inputplleindex,log10(inputlogEcut));
        sprintf(fgamname,"/Users/CEE/Desktop/NGC_253/pp_gam_spec_ss_finite_B_%.1f_plindex_%.3f_Ecut_%.1f.dat",(regBfield*1.e6),inputplleindex,log10(inputlogEcut));
        sprintf(fnuename,"/Users/CEE/Desktop/NGC_253/qe_nu_e_ss_finite_B_%.1f_plindex_%.3f_Ecut_%.1f.dat",(regBfield*1.e6),inputplleindex,log10(inputlogEcut));
        sprintf(fnumuname,"/Users/CEE/Desktop/NGC_253/qe_nu_mu_ss_finite_B_%.1f_plindex_%.3f_Ecut_%.1f.dat",(regBfield*1.e6),inputplleindex,log10(inputlogEcut));
    
        fne.open(fnename);
        fposi.open(fposiname);
        fgam.open(fgamname);
        fnue.open(fnuename);
        fnumu.open(fnumuname);
    
        cout << "Now doing photons and neutrinos.\n";
    
        logeps = logepsmin;
        while (logeps <= logepsmax) {
            eps = pow(10.0,logeps);
            qe_gamma = q_gamma(eps); // dN/dVdtdeps -- photons cm^-3 s^-1 eV^-1
            qe_nue_ss = q_nu_e(eps); // dN/dVdtdE -- neutrinos cm^-3 s^-1 eV^-1
            qe_numu_ss = q_nu_mu(eps); // dN/dVdtdE -- neutrinos cm^-3 s^-1 eV^-1
        
            //fgam << eps << " " << qe_gamma << endl;
        
            fgam << eps << "\t" << eps*eps*qe_gamma*volirregion*eV2ergs*1./(4.*pi*sourcedist*sourcedist) << endl; // ergs cm^-2 s^-1
            fnue << eps << "\t" << eps*eps*qe_nue_ss*volirregion*eV2ergs*1./(4.*pi*sourcedist*sourcedist) << endl; // ergs cm^-2 s^-1
            fnumu << eps << "\t" << eps*eps*qe_numu_ss*volirregion*eV2ergs*1./(4.*pi*sourcedist*sourcedist) << endl; // ergs cm^-2 s^-1
        
            cout << eps*(1.e-9) << " " << qe_gamma*eps*eps*(1.e-9) << endl; // GeV cm^-3 s^-1
        
            //cout << eps*(1.e-9) << " " << normalization*qe_gamma*1./(1.e-9) << endl; // photons cm^-3 s^-1 GeV^-1
        
            //fgam << eps << " " << normalization*eps*eps*qe_gamma << endl;
            //cout << logeps << " " << normalization*eps*eps*qe_gamma*eV2ergs << endl;
            //cout << logeps << " " << (normalization*(1.e9)/volirregion)*qe_gamma << endl;
        
            logeps += logepsbinsize;
        }
    
        fgam.close();
        fnue.close();
        fnumu.close();
    
        cout << "Now doing electrons.\n";
    
        loggamma_e = log_elec_gmin+logbinsize;
        while (loggamma_e <= log_elec_gmax) {
            E_elec = pow(10.0,loggamma_e)*elecmass; // eV
            qe_e_ss = q_elec(E_elec)*E_elec*log(10.0); // dN/dVdtd(log(gamma)) -- electrons cm^-3 s^-1 (log(gamma))^-1
            //qe_e_ss = normalization*q_e(E_elec)*0.5*E_elec*log(10.0); // dN/dVdtd(log(gamma)) -- electrons cm^-3 s^-1 (log(gamma))^-1
            //qe_p_ss = posi_spectrum_no_escape(loggamma_e);
            //qe_p_ss = qe_e_ss; // dN/dVdtd(log(gamma)) -- positrons cm^-3 s^-1 (log(gamma))^-1
            qe_p_ss = q_posi(E_elec)*E_elec*log(10.0); // dN/dVdtd(log(gamma)) -- electrons cm^-3 s^-1 (log(gamma))^-1
        
            cout << pow(10.0,loggamma_e) << " " << (qe_e_ss+qe_p_ss) << endl;
        
            //cout << pow(10.0,loggamma_e) << " " << elec_spectrum_noesc_pp(loggamma_e)*1./(loggamma_e*elecmass*(1.e-9)*log(10.0)) << endl; // dN/dVdtdE_e
        
            //cout << loggamma_e << " " << normalization*(qe_e_ss+qe_p_ss) << endl;
            //fne << loggamma_e << " " << normalization*(qe_e_ss+qe_p_ss) << endl;
        
            //cout << loggamma_e << " " << normalization*qe_e_ss << " " << normalization*qe_p_ss << endl;
        
            //fne << loggamma_e << " " << (normalization*1./volirregion)*qe_e_ss << endl;
            //fposi << loggamma_e << " " << (normalization*1./volirregion)*qe_p_ss << endl;
            
            fne << loggamma_e << "\t" << qe_e_ss << endl;
            fposi << loggamma_e <<"\t" << qe_p_ss << endl;
        
            //fne << loggamma_e << " " << normalization*qe_e_ss << endl;
            //fposi << loggamma_e << " " << normalization*qe_p_ss << endl;
        
            loggamma_e += logbinsize;
            //lines++;
        }
    
        //cout << lines << endl;
    
        fne.close();
        fposi.close();
        // inputplleindex loop ends here
    } 
    return 0;
}