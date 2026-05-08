#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <assert.h>

#include <stdbool.h>

#include "core_allvars.h"
#include "core_proto.h"

/* Diagnostic counters for the RK45 substep floor in integrate_spin_RK45().
 * Accumulated across the full run and printed by report_rk45_floor_stats(). */
static int    N_floor_triggers        = 0;   /* total floor invocations */
static int    N_calls_with_floor      = 0;   /* integrate_spin_RK45 calls that hit the floor >= once */
static double sum_stepfactor_at_floor = 0.0; /* for computing mean stepsize_factor at trigger */
static double min_stepfactor_at_floor = 1e30;
static double max_stepfactor_at_floor = 0.0;
static double sum_timefrac_at_first   = 0.0; /* time/dt_years at the FIRST floor hit per call */


//!! added declaration before implementation
double calcSpinupNT(double spin);

// Forward declaration for RK45 integrator
double integrate_spin_RK45(double *mass_ptr, double initial_spin, double initial_mdot_Edd,
                          double dt_years, int coherent, double madness,
                          int mad_superEddington, int mad_subEddington,
                          double allowedSpinError, double allowedFractionalMassError,
                          double minimumFractionalTimeResolution,
                          int *substeps_taken,
                          double *avg_jet_power, double *avg_magnetization,
                          double *avg_etaEM, double *avg_radiative_eff);



double r_isco(double spin) {
    double Z1 = 1.0 + pow(1.0 - pow(spin, 2), 1.0 / 3.0) * (pow(1.0 + spin, 1.0 / 3.0) + pow(1.0 - spin, 1.0 / 3.0));
    double Z2 = sqrt(3 * pow(spin, 2) + pow(Z1, 2));
    double r = 3 + Z2 - copysign(sqrt((3.0 - Z1) * (3.0 + Z1 + 2 * Z2)), spin);
    return r;
}

double radiativeEfficiencyThin(double spin) {
    return 1.0 - sqrt(1.0 - 2.0 / 3.0 / r_isco(spin));
}

double radiativeEfficiency(double spin, double mdot_Edd) {

    double hotAccretionSuppressionFactor = 100.0;
    double epsilon = radiativeEfficiencyThin(spin);
    
    if(mdot_Edd < fEdd_ADAF_threshold) {
        epsilon /= hotAccretionSuppressionFactor;
    } 
    else if(mdot_Edd > 5.0) {
        epsilon *= 5.0 / mdot_Edd * (1.0 + log(mdot_Edd / 5.0));
    }
    
    return epsilon;
}

double calcEddingtonLuminosity(double mass) {

    double L_Edd_cgs = 4.0 * M_PI * GRAVITY * mass * SOLAR_MASS * 
                       PROTONMASS * C / THOMSON_CROSS_SECTION;
    
    // Convert from erg/s to solar luminosities
    return L_Edd_cgs / SOLAR_LUM;
}

double calcEddingtonAccretionRateFiducial(double mass) {

    double fiducialRadiativeEfficiency = 0.1;
    
    // mdot_Edd = 4π G M m_p / (ε σ_T c) × yr
    return 4.0 * M_PI * GRAVITY * mass * SOLAR_MASS * PROTONMASS / 
           (fiducialRadiativeEfficiency * THOMSON_CROSS_SECTION * C) * 
           SEC_PER_YEAR / SOLAR_MASS;
}


double calcMaximumMagnetization(double spin) {
    return -20.2 * pow(spin, 3) - 14.9 * pow(spin, 2) + 34 * spin + 52.6;
}

double spinToHorizon(double a) {
    return 1.0 + sqrt(1.0 - pow(a, 2));
}

double Omega_H(double spin) {
    return fabs(spin) / (2 * spinToHorizon(spin));
}

double sHD_min(double spin) {
    return 0.86 - 2 * spin * 0.97;
}


double sHD(double spin, double r_Edd) {
    double parameters = 0.01971853;
    double xi = parameters * r_Edd;

    return (calcSpinupNT(spin) + sHD_min(spin) * xi) / (1 + xi);
}

double calcMagnetization(double spin, double r_Edd) {
    // Fitting parameters from RNC23
    const double param_0 = 1.29;
    const double param_1 = 1.88;

    double mdot_Edd = r_Edd * 0.1 / radiativeEfficiencyThin(spin);

    // Clamp inputs
    if(spin > 0.998) spin = 0.998;
    if(spin < -0.998) spin = -0.998;

    if(mdot_Edd > fEdd_ADAF_threshold) {
        // Thin disk regime: magnetization increases with accretion rate
        double evolutionFactor = pow(r_Edd / param_1, param_0);
        return calcMaximumMagnetization(spin) * evolutionFactor / (1.0 + evolutionFactor);
    } else {
        // Hot accretion regime: maximum magnetization (MAD state)
        return calcMaximumMagnetization(spin);
    }    

}

double jetEfficiencyFunction(double phi, double spin) {
    double kappa = 0.05;
    double horizonAngularVelocity = Omega_H(spin);
    return kappa / (4 * M_PI) * pow(phi, 2) * pow(horizonAngularVelocity, 2) * (1.0 + 1.38 * pow(horizonAngularVelocity, 2) - 9.2 * pow(horizonAngularVelocity, 4));
}

double etaEM(double spin, double r_Edd) {

    double efficiency = jetEfficiencyFunction(calcMagnetization(spin, r_Edd), spin);
    //if(efficiency > 0.99) efficiency = 0.99;  // Max 99% efficiency
    //if(efficiency < 0.0) efficiency = 0.0;    // Min 0%

    return efficiency;
}

double k_ratio(double spin) {
    if (spin <= 0) {
        return 0.23;
    } else {
        return fmin(0.35, 0.1 + 0.5 * spin);
    }
}


double calcJetPower(double mass, double spin, double mdot_Edd, 
                   double madness, int mad_superEddington, int mad_subEddington) {


    double massInflowRate = calcEddingtonAccretionRateFiducial(mass) * mdot_Edd;
    double r_Edd = mdot_Edd * radiativeEfficiencyThin(spin) / 0.1;
    int jets_enabled = (mad_superEddington && (mdot_Edd > fEdd_ADAF_threshold)) || 
                       (mad_subEddington && (mdot_Edd <= fEdd_ADAF_threshold));

    if(jets_enabled) {
        // Jet power = η_EM × Ṁ × c²
        double eta = etaEM(spin, r_Edd);
        
        // Convert mass inflow to grams/sec
        double massInflowRate_g_per_s = massInflowRate * SOLAR_MASS / SEC_PER_YEAR;
        
        // P_jet = η × Ṁ × c² in erg/s
        double jetPower_erg_per_s = eta * madness * massInflowRate_g_per_s * C * C;
        
        // Convert to Watts (1 W = 10^7 erg/s)
        return jetPower_erg_per_s / 1.0e7;
    } else {
        return 0.0;
    }

}

double calcMassEvolutionRate(double mass, double spin, double mdot_Edd, int coherent,
                            double madness, int mad_superEddington, int mad_subEddington) {

    double eta, eta_flip;
    double massInflowRate = calcEddingtonAccretionRateFiducial(mass) * mdot_Edd;
    double epsilon = radiativeEfficiency(spin, mdot_Edd);

    // Calculate jet efficiency
    int jets_enabled = (mad_superEddington && (mdot_Edd > fEdd_ADAF_threshold)) || 
                       (mad_subEddington && (mdot_Edd <= fEdd_ADAF_threshold));

    if(jets_enabled) {
        // eta = P_jet / (Ṁ × c²)
        double jetPower_watts = calcJetPower(mass, spin, mdot_Edd, madness, mad_superEddington, mad_subEddington);
        
        // Convert to erg/s
        double jetPower_erg_per_s = jetPower_watts * 1.0e7;
        
        // Mass-energy inflow rate: Ṁ × c² in erg/s
        double massInflowRate_g_per_s = massInflowRate * SOLAR_MASS / SEC_PER_YEAR;
        double massEnergyRate_erg_per_s = massInflowRate_g_per_s * C * C;
        
        eta = jetPower_erg_per_s / massEnergyRate_erg_per_s;
    } else {
        eta = 0.0;
    }

    // Loss factor: fraction of inflowing mass that stays on BH
    double lossFactor = 1.0 - epsilon - eta;

    if(!coherent) {
        // Incoherent: compute same for flipped orientation, then average
        double epsilon_flip = radiativeEfficiency(-spin, mdot_Edd);
        
        if(jets_enabled) {
            double jetPower_watts_flip = calcJetPower(mass, -spin, mdot_Edd, madness, mad_superEddington, mad_subEddington);
            double jetPower_erg_per_s_flip = jetPower_watts_flip * 1.0e7;
            double massInflowRate_g_per_s = massInflowRate * SOLAR_MASS / SEC_PER_YEAR;
            double massEnergyRate_erg_per_s = massInflowRate_g_per_s * C * C;
            eta_flip = jetPower_erg_per_s_flip / massEnergyRate_erg_per_s;
        } else{
            eta_flip = 0.0;
        }
        
        lossFactor += (1.0 - epsilon_flip - eta_flip);
        lossFactor /= 2.0;
    }

    // Return dM/dt in M☉/yr
    return massInflowRate * lossFactor;

}

double calcSpinupNT(double spin) {
    double a = fmax(-1, fmin(1, spin));
    double r = r_isco(spin);
    double E = (pow(r, 2) - 2 * r + a * sqrt(r)) / (r * sqrt(pow(r, 2) - 3 * r + 2 * a * sqrt(r)));
    double l = sqrt(r) * (pow(r, 2) - 2 * a * sqrt(r) + pow(a, 2)) / (r * sqrt(pow(r, 2) - 3 * r + 2 * a * sqrt(r)));
    return l - 2 * a * E;
}

double calcSpinupMAD(double spin) {
    double fitCoefficients[6] = {0.44763054, -12.53175763, -7.79639428, 9.43649497, 5.70845519, -4.02630245};
    double result = 0;
    for (int i = 0; i < 6; ++i) {
        result += fitCoefficients[i] * pow(spin, i);
    }

    return result;
}


double calcSpinup(double spin, double mdot_Edd, double madness, int mad_superEddington, int mad_subEddington) {

    double r_Edd = mdot_Edd * radiativeEfficiencyThin(spin) / 0.1;
    double s;
    
    if(mad_subEddington && (mdot_Edd < fEdd_ADAF_threshold)) {
        // MAD ADAF regime (hot accretion)
        s = calcSpinupMAD(spin) * madness + calcSpinupNT(spin) * (1.0 - madness);
    }
    else if(mad_superEddington && (mdot_Edd >= fEdd_ADAF_threshold)) {
        // Thin disk regime (super-Eddington relative to ADAF threshold)
        
        // Start with hydrodynamic spinup
        s = sHD(spin, r_Edd);
        
        // Add electromagnetic spinup if spin is non-zero
        if(spin != 0.0) {
            double eta = etaEM(spin, r_Edd);
            double k_rat = k_ratio(spin);
            double omega_H = Omega_H(spin);
            double spin_sign = (spin > 0.0) ? 1.0 : -1.0;
            
            s += -eta * (1.0 / k_rat / omega_H - 2.0 * spin) * spin_sign;
        }
        
        // Blend with Novikov-Thorne if madness < 1.0
        if(madness < 1.0) {
            s = s * madness + calcSpinupNT(spin) * (1.0 - madness);
        }
    }
    else {
        // Default to Novikov-Thorne (no MAD physics)
        s = calcSpinupNT(spin);
    }
    
    return s;
}

double calcSpinEvolutionRate(double mass, double spin, double mdot_Edd, int coherent,
                            double madness, int mad_superEddington, int mad_subEddington) {

    double s_coherent = calcSpinup(spin, mdot_Edd, madness, mad_superEddington, mad_subEddington);
    double mdot_Edd_fiducial = calcEddingtonAccretionRateFiducial(mass);
    double coherentSpinup = s_coherent * mdot_Edd * mdot_Edd_fiducial / mass;

    if(coherent) {
        // Protect against extreme spin evolution rates
        // if(coherentSpinup != coherentSpinup) coherentSpinup = 0.0;  // NaN
        // if(coherentSpinup > 100.0) coherentSpinup = 100.0;   // Cap at 100/yr
        // if(coherentSpinup < -100.0) coherentSpinup = -100.0;
        return coherentSpinup;
    } else {
        // Assume that spinflips happen all the time, on timescales short compared 
        // to the evolutionary timescale. Careful--there's a sign flip convention required.
        double s_incoherent = calcSpinup(-spin, mdot_Edd, madness, mad_superEddington, mad_subEddington);
        double incoherentSpinup = -s_incoherent * mdot_Edd * mdot_Edd_fiducial / mass;
        double spin_av = (coherentSpinup + incoherentSpinup) / 2.0;
        // Protect against extreme spin evolution rates
        // if(spin_av != spin_av) spin_av = 0.0;  // NaN
        // if(spin_av > 100.0) spin_av = 100.0;   // Cap at 100/yr
        // if(spin_av < -100.0) spin_av = -100.0;

        return spin_av;
    }

}


//!! added stadard coefficients definition .. it makes no sense to put them in struct at each instantiation of the struct
const struct {
    double _A[6];
    double _B[6][5]; 
    double _C[6] ;
    double _CH[6];
    double _CT[6];  
}RKFCoefficients = {
    {0.0, 2.0 / 9.0, 1.0 / 3.0, 3.0 / 4.0, 1.0, 5.0 / 6.0},//_A
    {{0.0, 0.0, 0.0, 0.0, 0.0}, 
     {2.0 / 9.0, 0.0, 0.0, 0.0, 0.0}, 
     {1.0 / 12.0, 1.0 / 4.0, 0.0, 0.0, 0.0}, 
     {69.0 / 128.0, -243.0 / 128.0, 135.0 / 64.0, 0.0, 0.0}, 
     {-17.0 / 12.0, 27.0 / 4.0, -27.0 / 5.0, 16.0 / 15.0, 0.0}, 
     {65.0 / 432.0, -5.0 / 16.0, 13.0 / 16.0, 4.0 / 27.0, 5.0 / 144.0}},//_B
    {1.0 / 9.0, 0.0, 9.0 / 20.0, 16.0 / 45.0, 1.0 / 12.0},//_C
    {47.0 / 450.0, 0.0, 12.0 / 25.0, 32.0 / 225.0, 1.0 / 30.0, 6.0 / 25.0}, //_CH
    {-1.0 / 150.0, 0.0, 3.0 / 100.0, -16.0 / 75.0, -1.0 / 20.0, 6.0 / 25.0} //_CT   
};



double integrate_spin_RK45(double *mass_ptr, double initial_spin, double initial_mdot_Edd,
                          double dt_years, int coherent, double madness,
                          int mad_superEddington, int mad_subEddington,
                          double allowedSpinError, double allowedFractionalMassError,
                          double minimumFractionalTimeResolution,
                          int *substeps_taken,
                          double *avg_jet_power, double *avg_magnetization,
                          double *avg_etaEM, double *avg_radiative_eff)
{
    double mass = *mass_ptr;
    double spin = initial_spin;
    double time = 0.0;
    double substep = dt_years * 0.01;  // Initial substep (1% of full timestep)

    double acc_jet_power    = 0.0;
    double acc_magnetization = 0.0;
    double acc_etaEM        = 0.0;
    double acc_radiative_eff = 0.0;
    double acc_time         = 0.0;

    
    const int MAX_SUBSTEPS = 500;  // Hard limit as requested
    int substep_count = 0;
    int this_call_hit_floor = 0;  /* set on the first floor trigger in this call */

    // Compute actual mass inflow rate (Msun/yr) — constant over timestep
    double mdot_actual = initial_mdot_Edd * calcEddingtonAccretionRateFiducial(*mass_ptr); 
    // Integrate until we reach dt_years
    while(time < dt_years) {

        if(substep_count >= MAX_SUBSTEPS) {
            break;
        }        

        // Don't overshoot the target time
        if(time + substep > dt_years) {
            substep = dt_years - time;
        }      
        //

        int step_accepted = 0;
        double proposed_spin, proposed_mass;       

        while(!step_accepted) {
            // RK45 substep calculations
            double k_spin[6], k_mass[6];
            
            for(int k = 0; k < 6; k++) {
                double mass_pred = mass;
                double spin_pred = spin;
                
                for(int l = 0; l < k; l++) {
                    mass_pred += k_mass[l] * RKFCoefficients._B[k][l];
                    spin_pred += k_spin[l] * RKFCoefficients._B[k][l];
                }
                
                // Enforce spin bounds during prediction
                if(spin_pred > 0.998) spin_pred = 0.998;
                if(spin_pred < -0.998) spin_pred = -0.998;
                if(mass_pred < 1e-10) mass_pred = 1e-10;

                // check with Angelo
                // Now recompute mdot_Edd at predicted mass
                double mdot_Edd_pred = mdot_actual / calcEddingtonAccretionRateFiducial(mass_pred);
                // check with Angelo

                // Calculate da/dt at this prediction point
                double dM_dt = calcMassEvolutionRate(mass_pred, spin_pred, mdot_Edd_pred, coherent,
                                                     madness, mad_superEddington, mad_subEddington);
                double da_dt = calcSpinEvolutionRate(mass_pred, spin_pred, mdot_Edd_pred, coherent,
                                                     madness, mad_superEddington, mad_subEddington);

                k_mass[k] = substep * dM_dt;
                k_spin[k] = substep * da_dt;
            }
            
            // Proposed new spin (5th order)
            proposed_mass = mass;
            proposed_spin = spin;
            for(int k = 0; k < 6; k++) {
                proposed_mass += RKFCoefficients._CH[k] * k_mass[k];
                proposed_spin += RKFCoefficients._CH[k] * k_spin[k];
            }
            
            // Enforce bounds
            if(proposed_spin > 0.998) proposed_spin = 0.998;
            if(proposed_spin < -0.998) proposed_spin = -0.998;
            if(proposed_mass < 1e-10) proposed_mass = 1e-10;
            
            // Estimate truncation error
            double truncation_error_mass = 0.0;
            double truncation_error_spin = 0.0;
            for(int k = 0; k < 6; k++) {
                truncation_error_mass += RKFCoefficients._CT[k] * k_mass[k];
                truncation_error_spin += RKFCoefficients._CT[k] * k_spin[k];
            }
            truncation_error_mass = fabs(truncation_error_mass);
            truncation_error_spin = fabs(truncation_error_spin);

            double stepsize_factor_mass, stepsize_factor_spin;

            stepsize_factor_mass = (truncation_error_mass > 0) ? 0.9 * pow(allowedFractionalMassError / (truncation_error_mass / mass), 0.2) : 1e30;
            stepsize_factor_spin = (truncation_error_spin > 0) ? 0.9 * pow(allowedSpinError / truncation_error_spin, 0.2) : 1e30;

            //stepsize_factor_mass = 0.9 * pow(allowedFractionalMassError / (truncation_error_mass / mass), 0.2);
            //stepsize_factor_spin = 0.9 * pow(allowedSpinError / truncation_error_spin, 0.2);

            // double stepsize_factor_spin = 0.9 * pow(allowedSpinError / (truncation_error + 1e-30), 0.2);
            
            // Accept step if within tolerance (matches Python line 368)
            if(stepsize_factor_mass >= 1.0 && stepsize_factor_spin >= 1.0) {
                step_accepted = 1;
                mass = proposed_mass;
                spin = proposed_spin;
                time += substep;
                substep_count++;

                // accumulate time-weighted averages 
                double mdot_Edd_step = mdot_actual / calcEddingtonAccretionRateFiducial(mass);
                double r_Edd_step    = mdot_Edd_step * radiativeEfficiencyThin(spin) / 0.1;
                acc_jet_power     += calcJetPower(mass, spin, mdot_Edd_step, madness,
                                                  mad_superEddington, mad_subEddington) * substep;
                acc_magnetization += calcMagnetization(spin, r_Edd_step) * substep;
                acc_etaEM         += etaEM(spin, r_Edd_step)             * substep;
                acc_radiative_eff += radiativeEfficiencyThin(spin)        * substep;
                acc_time          += substep;

            }
           
            // Adjust timestep (matches Python line 372)
            // Limit increase to factor of 100
            double stepsize_factor = fmin(stepsize_factor_mass, stepsize_factor_spin);
            substep *= fmin(stepsize_factor, 100.0);

            // Ensure minimum resolution
            if( (minimumFractionalTimeResolution > 0) && (substep > dt_years * minimumFractionalTimeResolution) ) {
                substep = dt_years * minimumFractionalTimeResolution;
            }

            // Floor: if substep has shrunk below machine-precision threshold the inner loop
            // would stall indefinitely (k values -> 0, error -> 0, step "accepted" with zero
            // time advance). Force-accept the current proposed state at the floor size.
            // This can occur for very stiff ODEs encountered in massive TNG100 halos.
            if (!step_accepted && substep < dt_years * 1e-10) {
                /* Diagnostics: record how stiff the ODE was and how far we had integrated. */
                N_floor_triggers++;
                if (!this_call_hit_floor) {
                    this_call_hit_floor = 1;
                    N_calls_with_floor++;
                    sum_timefrac_at_first += time / dt_years;
                }
                if (stepsize_factor < min_stepfactor_at_floor) min_stepfactor_at_floor = stepsize_factor;
                if (stepsize_factor > max_stepfactor_at_floor) max_stepfactor_at_floor = stepsize_factor;
                sum_stepfactor_at_floor += stepsize_factor;

                substep    = dt_years * 1e-10;
                step_accepted = 1;
                mass = proposed_mass;
                spin = proposed_spin;
                time += substep;
                substep_count++;
                double mdot_Edd_step = mdot_actual / calcEddingtonAccretionRateFiducial(mass);
                double r_Edd_step    = mdot_Edd_step * radiativeEfficiencyThin(spin) / 0.1;
                acc_jet_power     += calcJetPower(mass, spin, mdot_Edd_step, madness,
                                                  mad_superEddington, mad_subEddington) * substep;
                acc_magnetization += calcMagnetization(spin, r_Edd_step) * substep;
                acc_etaEM         += etaEM(spin, r_Edd_step)             * substep;
                acc_radiative_eff += radiativeEfficiencyThin(spin)        * substep;
                acc_time          += substep;
            }
        }
    }

    // write time-averaged outputs 
    if(acc_time > 0.0) {
        if(avg_jet_power)     *avg_jet_power     = acc_jet_power     / acc_time;
        if(avg_magnetization) *avg_magnetization = acc_magnetization / acc_time;
        if(avg_etaEM)         *avg_etaEM         = acc_etaEM         / acc_time;
        if(avg_radiative_eff) *avg_radiative_eff = acc_radiative_eff / acc_time;
    } else {
        // no steps completed — fallback to final state
        double mdot_f = mdot_actual / calcEddingtonAccretionRateFiducial(mass);
        double r_Edd_f = mdot_f * radiativeEfficiencyThin(spin) / 0.1;
        if(avg_jet_power)     *avg_jet_power     = calcJetPower(mass, spin, mdot_f, madness,
                                                   mad_superEddington, mad_subEddington);
        if(avg_magnetization) *avg_magnetization = calcMagnetization(spin, r_Edd_f);
        if(avg_etaEM)         *avg_etaEM         = etaEM(spin, r_Edd_f);
        if(avg_radiative_eff) *avg_radiative_eff = radiativeEfficiencyThin(spin);
    }

    // Return number of substeps taken
    if(substeps_taken != NULL) *substeps_taken = substep_count;
    *mass_ptr = mass;  // Update final mass
    return spin;
}



// Simple wrapper function to update spin during accretion
void update_black_hole_spin(int p, double accreted_mass, double dt, int coherent,
                           double madness, int mad_superEddington, int mad_subEddington) {

    double mass_msun, spin, mdot_edd_initial, mdot_actual, mdot_Edd_initial, new_spin, new_mass_msun;
    double mdot_Edd_final, magnetization, jet_power, jet_efficiency, radiative_eff;
    double avg_jet_power = 0.0, avg_magnetization = 0.0;
    double avg_etaEM = 0.0, avg_radiative_eff = 0.0;


    mass_msun = Gal[p].BlackHoleMass * 1e10 / Hubble_h;
    spin = Gal[p].BlackHoleSpin;
    
    // safety check
    // Clamp spin BEFORE any calculations
    if(spin > 0.998) spin = 0.998;
    if(spin < -0.998) spin = -0.998;
    if(spin != spin) spin = SpinSeedValue;  // NaN check
    
    // Calculate INITIAL Eddington ratio (for RK45 integration)
    mdot_edd_initial = calcEddingtonAccretionRateFiducial(mass_msun);
    mdot_actual = (accreted_mass * 1e10 / Hubble_h) / (dt * UnitTime_in_s / SEC_PER_YEAR);
    mdot_Edd_initial = mdot_actual / mdot_edd_initial;
    
    if(mdot_Edd_initial > 1000.0) mdot_Edd_initial = 1000.0;


    // RK45 parameters
    double allowedSpinError = 1e-5;
    double allowedFractionalMassError = 1e-3;
    double minimumFractionalTimeResolution = 1.0;
    double dt_years = dt * UnitTime_in_s / SEC_PER_YEAR; 

    // Track substeps
    int substeps_taken = 0;

    // Integrate spin using RK45
    new_mass_msun = mass_msun;
    new_spin = integrate_spin_RK45(&new_mass_msun, spin, mdot_Edd_initial, dt_years, coherent,
                               madness, mad_superEddington, mad_subEddington,
                               allowedSpinError, allowedFractionalMassError,
                               minimumFractionalTimeResolution, &substeps_taken,
                               &avg_jet_power, &avg_magnetization,
                               &avg_etaEM, &avg_radiative_eff);
    

    // For incoherent: randomly choose prograde or retrograde (matches Python line 381-382)
    if(!coherent) {
        new_spin = (drand48() < 0.5) ? fabs(new_spin) : -fabs(new_spin);
    }
    
    double DeltaBlackHoleMass = (new_mass_msun * Hubble_h / 1e10) - Gal[p].BlackHoleMass;
    Gal[p].DeltaBlackHoleMass = DeltaBlackHoleMass;  // net increment (after losses) from this sub-step's accretion
    Gal[p].BlackHoleMass = new_mass_msun * Hubble_h / 1e10;
    Gal[p].BlackHoleSpin = new_spin;
 
    // Calculate final properties using NEW mdot_Edd
    magnetization  = avg_magnetization;
    jet_power      = avg_jet_power;
    jet_efficiency = avg_etaEM;
    radiative_eff  = avg_radiative_eff;

    // Store final values
    Gal[p].Eddrate = mdot_Edd_initial;  // Use NEW mdot_Edd
    Gal[p].BH_Magnetization = magnetization;
    Gal[p].BH_JetPower = jet_power;
    Gal[p].BH_JetEfficiency = jet_efficiency;
    Gal[p].BH_thindiskRadiativeEfficiency = radiative_eff;

    // DEBUG: print spin and accretion rates if BH mass more than doubled this sub-step, then crash
    if(new_mass_msun / mass_msun > 2.0)
    {
        printf("DEBUG update_black_hole_spin: spin_initial=%.6f  spin_final=%.6f  "
               "mdot_Edd_initial=%.6e  mdot_actual=%.6e (Msun/yr)\n",
               spin, new_spin, mdot_Edd_initial, mdot_actual);
        fflush(stdout);
        abort();
    }
}

void update_black_hole_spin_coherent(int p, double accreted_mass, double dt) {
    update_black_hole_spin(p, accreted_mass, dt, 1, 1.0, 1, 1);  // coherent, madness=1.0, both regimes enabled
}

void update_black_hole_spin_incoherent(int p, double accreted_mass, double dt) {
    update_black_hole_spin(p, accreted_mass, dt, 0, 1.0, 1, 1);  // incoherent, madness=1.0, both regimes enabled
}


// evolve_black_hole: the single unified BH evolution call, made once per sub-step per galaxy
// in core_build_model.c.
//
// It reads BH_incoherent_reservoir (hot-gas / radio-mode supply accumulated this sub-step) and
// BH_coherent_reservoir (cold-disc / disk-instability supply accumulated this sub-step), resets
// both to zero, then performs ONE call to update_black_hole_spin() using the combined total.
//
// Whether the spin update is coherent or incoherent is determined by whichever reservoir is
// larger -- if the disk-instability supply dominates, accretion is treated as coherent (aligned
// with the disc angular momentum); if the hot-gas supply dominates, it is incoherent (random
// orientation flips on short timescales).
//
// After the spin update the net mass increment (DeltaBlackHoleMass, after radiative + jet losses)
// is split back into radio and instability fractions proportionally by supply mass, so that
// RadioBlackHoleMass, InstaBlackHoleMass, etc. are tracked correctly.
//
// Returns: radio_net -- the net BH mass that came from the radio-mode (hot-gas) supply.
//          This is passed on to bh_feedback() so it can compute the AGN jet heating suppression.
//
// DeltaBlackHoleMass is set to the instability-driven fraction of the net increment so that
// core_build_model.c can do:  InstaBlackHoleMass += DeltaBlackHoleMass  there.
double evolve_black_hole(int p, double dt)
{
    double radio_supply  = Gal[p].BH_incoherent_reservoir;
    double instab_supply = Gal[p].BH_coherent_reservoir;
    double total_supply  = radio_supply + instab_supply;

    // Zero the reservoirs now that we've read them; they will refill during the next sub-step.
    Gal[p].BH_incoherent_reservoir = 0.0;
    Gal[p].BH_coherent_reservoir   = 0.0;

    if(total_supply <= 0.0)
    {
        // Nothing to accrete this sub-step.  Zero the staging field so that
        // InstaBlackHoleMass is not inflated by a stale value from a previous sub-step.
        Gal[p].DeltaBlackHoleMass = 0.0;
        return 0.0;
    }

    // Fraction of the total supply that came from the radio-mode (hot-gas) channel.
    // Used below to split the net BH mass increment back into radio vs. instability bookkeeping.
    double radio_fraction = radio_supply / total_supply;
    double radio_net = 0.0;

    if(SpinEvolutionOn > 0)
    {
        // Single unified spin update using the combined supply from both channels.
        // Coherent (disc-aligned) accretion applies when the instability supply dominates;
        // incoherent (random-orientation) accretion applies when the hot-gas supply dominates.
        int coherent = (instab_supply >= radio_supply) ? 1 : 0;
        update_black_hole_spin(p, total_supply, dt, coherent, 1.0, 1, 1);

        // update_black_hole_spin() has updated BlackHoleMass and set DeltaBlackHoleMass to
        // the total net increment (after radiative + jet losses) from both channels combined.
        double total_net = Gal[p].DeltaBlackHoleMass;

        // Apportion the net increment between the two supply channels proportionally.
        // BlackHoleMass itself is already correct; we only need to update the bookkeeping fields.
        radio_net         = total_net * radio_fraction;
        double instab_net = total_net - radio_net;

        Gal[p].RadioBlackHoleMass += radio_net;
        Gal[p].RadioBHaccreted     = radio_net / dt;
        if(radio_net / dt > Gal[p].MaxRadioModeAccretionRate)
            Gal[p].MaxRadioModeAccretionRate = radio_net / dt;

        // Leave DeltaBlackHoleMass holding only the instability fraction so that
        // core_build_model.c can do:  InstaBlackHoleMass += DeltaBlackHoleMass
        Gal[p].DeltaBlackHoleMass = instab_net;

        // Accumulate mass lost to radiation and jets this sub-step.
        // total_net already accounts for both radiative and jet losses inside the RK45 integrator.
        Gal[p].BH_CumulativeMassLoss += (total_supply - total_net);
    }
    else
    {
        // SpinEvolutionOn == 0: simple fixed-efficiency accretion, no spin physics.

        // Radio-mode: apply the global RadiativeEfficiency to get the net BH mass gain.
        // (The accreted supply was already removed from HotGas in cooling_recipe_part_a.)
        radio_net = (1.0 - RadiativeEfficiency) * radio_supply;
        Gal[p].BlackHoleMass   += radio_net;
        Gal[p].BondiBHaccreted  = radio_net / dt;  // Bondi-style radio-mode accretion rate
        Gal[p].RadioBlackHoleMass += radio_net;
        Gal[p].RadioBHaccreted  = radio_net / dt;
        if(radio_net / dt > Gal[p].MaxRadioModeAccretionRate)
            Gal[p].MaxRadioModeAccretionRate = radio_net / dt;

        // Disk instability: in the old code the gas went directly to BlackHoleMass with no
        // radiative efficiency applied, so we preserve that behaviour here.
        Gal[p].BlackHoleMass += instab_supply;

        // Record only the instability-driven increment in DeltaBlackHoleMass so that
        // InstaBlackHoleMass is updated correctly in core_build_model.c.
        // (The old code erroneously used the radio-mode increment here; this fixes that bug.)
        Gal[p].DeltaBlackHoleMass = instab_supply;
    }

    assert(Gal[p].BlackHoleMass >= 0.0);
    return radio_net;
}


// bh_feedback: unified AGN feedback, applied once per sub-step per galaxy in core_build_model.c
// immediately after evolve_black_hole().
//
// Two feedback channels are handled here:
//
//   1. Radio-mode (jet) heating -- suppresses cooling gas by heating it back toward the hot-halo
//      energy state.  For SpinEvolutionOn==0 this uses the classic AGNcoeff formula
//      (8.98755e10 * epsilon / specific_energy_change); for SpinEvolutionOn>0 it uses the
//      time-averaged BH_JetPower set by the spin update that was just completed inside
//      evolve_black_hole().  Both cases compute AGNheating first, then cap it at coolingGas,
//      and set AGNfountain to the *uncapped* value so the fountain can exceed the cooling rate.
//
//   2. Quasar-mode winds -- driven by instab_supply (the disk-instability BH supply).  Calling
//      this here -- after evolve_black_hole() has updated BH_thindiskRadiativeEfficiency via the
//      spin integrator -- fixes a bug in the old code where quasar_mode_wind() was invoked in
//      check_disk_instability() *before* update_black_hole_spin_coherent(), meaning it used the
//      radio-mode spin efficiency rather than the disk-instability one.
//
// The fountain code is preserved verbatim from the old do_AGN_heating().
//
// Returns: suppressed coolingGas (= coolingGas_raw minus AGN heating suppression).
double bh_feedback(int p, double radio_supply, double instab_supply, double radio_net,
                   double coolingGas_raw, double specific_energy_change, double dt)
{
    double coolingGas  = coolingGas_raw;
    double AGNheating  = 0.0;
    double AGNfountain = 0.0;

    // Radio-mode AGN heating: only meaningful when there was a hot-gas supply and a non-zero
    // specific energy difference (i.e. the hot gas is energetically above the cold disc state).
    // specific_energy_change == 0 implies either no hot gas or stable gas, so we skip.
    if(AGNrecipeOn > 0 && (radio_supply + instab_supply) > 0.0 && specific_energy_change > 0.0)
    {
        if(SpinEvolutionOn == 0)
        {
            // coefficient to heat the cooling gas back to the energy state of the hot halo
            // 8.98755e10 = c^2 in (km/s)^2
            double AGNcoeff = 8.98755e10 * RadiativeEfficiency / specific_energy_change;

            // cooling mass that can be suppressed from AGN heating
            AGNheating  = AGNcoeff * radio_supply;
            AGNfountain = AGNheating; // set BEFORE the cap below so fountain can exceed cooling rate
        }
        else // SpinEvolutionOn > 0
        {
            // BH_JetPower was set by update_black_hole_spin() inside evolve_black_hole(),
            // reflecting the time-averaged jet output from the RK45 spin integration.
            double dt_sec            = dt * UnitTime_in_s;
            double jet_energy        = Gal[p].BH_JetPower * 1e7 * dt_sec;    // erg  (1 W = 1e7 erg/s)
            double eta_floor= 0.1;
            double maccreted_g = (radio_supply + instab_supply) * 1e10 / Hubble_h * SOLAR_MASS; // g
            double jet_energy_floor = eta_floor * maccreted_g * C * C;        // erg
            if(jet_energy < jet_energy_floor) jet_energy = jet_energy_floor;
            double specific_energy_cgs = specific_energy_change * 1e10;       // (km/s)^2 --> erg/g

            // cooling mass that can be suppressed from AGN heating
            AGNheating  = EpsilonR * jet_energy / specific_energy_cgs / (1e10 * SOLAR_MASS / Hubble_h);
            AGNfountain = AGNheating; // set BEFORE the cap below so fountain can exceed cooling rate
        }

        // limit heating to cooling rate.
        // Making a conscious decision to no longer update the AGN accretion rate in proportion
        // here.  In effect, if this if-statement is triggered, the energy coupling for the radio
        // mode temporarily decreases.
        if(AGNheating > coolingGas)
            AGNheating = coolingGas;

        // Recalculate the actual radio-mode accretion rate to potentially update the historical
        // maximum.  We use radio_net (the post-loss net increment) consistent with how
        // evolve_black_hole() tracks MaxRadioModeAccretionRate.
        double AGNrate = radio_net / dt;
        if(AGNrate > Gal[p].MaxRadioModeAccretionRate)
            Gal[p].MaxRadioModeAccretionRate = AGNrate;

        Gal[p].AGNheating  = AGNheating;
        Gal[p].AGNfountain = AGNfountain;
        Gal[p].Heating    += AGNheating * specific_energy_change; // energy from the AGN pumped into keeping the hot gas hot

        // Fountain code: gas kicked up by the AGN jet is moved from HotGas into FountainGas.
        // AGNfountain may exceed AGNheating (and hence coolingGas) because it is the uncapped
        // estimate; the fountain can mobilise more mass than the effective cooling suppression.
        double metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
        if(AGNfountain > 0.0)
        {
            if(AGNfountain >= Gal[p].HotGas)
            {
                Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + Gal[p].HotGas * 0.1 / sqrt(Hubble_sqr_z(Halo[Gal[p].HaloNr].SnapNum))) / (Gal[p].FountainGas + Gal[p].HotGas);
                Gal[p].FountainGas      += Gal[p].HotGas;
                Gal[p].MetalsFountainGas += Gal[p].MetalsHotGas;
                Gal[p].HotGas       = 0.0;
                Gal[p].MetalsHotGas = 0.0;
            }
            else
            {
                Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + AGNfountain * 0.1 / sqrt(Hubble_sqr_z(Halo[Gal[p].HaloNr].SnapNum))) / (Gal[p].FountainGas + AGNfountain);
                Gal[p].FountainGas      += AGNfountain;
                Gal[p].MetalsFountainGas += (metallicity * AGNfountain);
                Gal[p].HotGas       -= AGNfountain;
                Gal[p].MetalsHotGas -= (metallicity * AGNfountain);
                if(Gal[p].MetalsHotGas < 0.0) Gal[p].MetalsHotGas = 0.0;
 		if(Gal[p].MetalsFountainGas < 0.0) Gal[p].MetalsFountainGas = 0.0;
            }
        }

        coolingGas -= AGNheating;
    }

    // Quasar-mode wind: driven by the disk-instability BH supply (instab_supply).
    // Calling this here -- after evolve_black_hole() has updated BH_thindiskRadiativeEfficiency
    // via the spin integrator -- ensures the correct radiative efficiency is used.
    // In the old code quasar_mode_wind() was called in check_disk_instability() *before*
    // update_black_hole_spin_coherent(), which was a bug.
    if((radio_supply+instab_supply) > 0.0 && AGNrecipeOn > 0)
        quasar_mode_wind(p, (float)(radio_supply+instab_supply));

    // Guard against floating-point rounding driving coolingGas slightly negative.
    if(coolingGas < 0.0) coolingGas = 0.0;

    return coolingGas;
}


double calculate_merger_spin_aligned(double M1, double a1, double M2, double a2,
                                     double theta1, double theta2, 
                                     double phi1, double phi2)
{
    // Final BH spin from aligned merger
    double q, nu, a_final, l_orb;
    double s4, s5, t0, t2, t3;
    double cosAlpha, cosBeta, cosGamma;
    double underRadical, J_total;
    
    // Ensure M1 >= M2 (swap if needed)
    if(M2 > M1) {
        double temp_M = M1, temp_a = a1, temp_theta = theta1, temp_phi = phi1;
        M1 = M2; a1 = a2; theta1 = theta2; phi1 = phi2;
        M2 = temp_M; a2 = temp_a; theta2 = temp_theta; phi2 = temp_phi;
    }
    
    // Handle edge cases
    if(M1 < 1e-10) return 0.0;  // Both essentially zero
    if(M2 < 1e-10) return a1;              // Only M1 has mass
    
    // Mass ratio q = M2/M1 <= 1
    q = M2 / M1;
    
    // Symmetric mass ratio nu = q/(1+q)^2
    nu = q / ((1.0 + q) * (1.0 + q));
    
    // Fitting coefficients from Rezzolla et al. (2008) Table 1
    // These come from fitting to numerical relativity simulations
    s4 = -0.129;
    s5 = -0.384;
    t0 = -2.686;
    t2 = -3.454;
    t3 = 2.353;
    
    // Compute angles between spins
    cosBeta = cos(theta1);   // Angle between a1 and orbital angular momentum
    cosGamma = cos(theta2);  // Angle between a2 and orbital angular momentum
    
    // Angle between a1 and a2 (dot product of spin vectors)
    if(a1 > 1e-10 && a2 > 1e-10) {
        // a_vec_1 · a_vec_2 = a1*a2*cos(angle between them)
        cosAlpha = sin(theta1)*sin(theta2)*cos(phi1 - phi2) + cos(theta1)*cos(theta2);
    } else {
        cosAlpha = 1.0;
    }

    // Specific orbital angular momentum at merger
    // This captures the contribution from the orbital motion
    l_orb = s4 / ((1.0 + q*q) * (1.0 + q*q)) * (a1*a1 + a2*a2*q*q*q*q + 2.0*a1*a2*q*q*cosAlpha)
            + (s5*nu + t0 + 2.0) / (1.0 + q*q) * (a1*cosBeta + a2*q*q*cosGamma)
            + 2.0*sqrt(3.0) + t2*nu + t3*nu*nu;
 
    // Total angular momentum in orbital direction (z-component)
    J_total = a1*cosBeta + a2*q*q*cosGamma + l_orb*q;

    // Final spin formula (Rezzolla et al. 2008, Eq. 7)
    a_final = sqrt(a1*a1 + a2*a2*q*q*q*q + 2.0*a2*a1*q*q*cosAlpha
                   + 2.0*(a1*cosBeta + a2*q*q*cosGamma)*l_orb*q
                   + l_orb*l_orb*q*q) / ((1.0 + q) * (1.0 + q));
   
    // Apply sign based on total angular momentum direction
    if(J_total < 0) {
        a_final = -a_final;  // Retrograde merger
    }    

    // For near-equal mass mergers (q > 0.05), randomly flip sign with 50% probability
    if(q > 0.05 && drand48() > 0.5) {
        a_final = -a_final;
    }

    // Enforce physical limits: -0.998 < a < 0.998
    // (Using 0.998 instead of 1.0 to match your other spin limits)
    if(a_final > 0.998) a_final = 0.998;
    if(a_final < -0.998) a_final = -0.998;

    return a_final;
}


void report_rk45_floor_stats(void)
{
    printf("\n--- RK45 substep floor diagnostics ---\n");
    if (N_floor_triggers == 0) {
        printf("  No floor triggers recorded.\n");
        printf("--------------------------------------\n\n");
        return;
    }
    printf("  Total floor triggers:                    %d\n", N_floor_triggers);
    printf("  integrate_spin_RK45 calls that hit floor: %d\n", N_calls_with_floor);
    printf("  Mean stepsize_factor at trigger:         %.4e  (1.0 = converged)\n",
           sum_stepfactor_at_floor / N_floor_triggers);
    printf("  Min stepsize_factor at trigger:          %.4e  (most stiff)\n",  min_stepfactor_at_floor);
    printf("  Max stepsize_factor at trigger:          %.4e  (closest to converged)\n", max_stepfactor_at_floor);
    printf("  Mean time fraction covered at first hit: %.4f  (0=stalled immediately, 1=near end)\n",
           sum_timefrac_at_first / N_calls_with_floor);
    printf("--------------------------------------\n\n");
}
