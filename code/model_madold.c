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


//!! added declaration before implementation
double calcSpinupNT(double spin);

// Forward declaration for RK45 integrator
double integrate_spin_RK45(double *mass_ptr, double initial_spin, double initial_mdot_Edd,
                          double dt_years, int coherent, double madness,
                          int mad_superEddington, int mad_subEddington,
                          double allowedSpinError, double allowedFractionalMassError,
                          double minimumFractionalTimeResolution,
                          int *substeps_taken);


double r_isco(double spin) {
    double Z1 = 1.0 + pow(1.0 - pow(spin, 2), 1.0 / 3.0) * (pow(1.0 + spin, 1.0 / 3.0) + pow(1.0 - spin, 1.0 / 3.0));
    double Z2 = sqrt(3 * pow(spin, 2) + pow(Z1, 2));
    double r = 3 + Z2 - copysign(sqrt((3.0 - Z1) * (3.0 + Z1 + 2 * Z2)), spin);
    return r;
}

double radiativeEfficiencyThin(double spin) {
    return 1.0 - sqrt(1.0 - 2.0 / 3.0 / r_isco(spin));
}

double radiativeEfficiency(double spin, double eddingtonRatio) {

    double mdot_Edd = eddingtonRatio * 0.1 / radiativeEfficiencyThin(spin);
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


double sHD(double spin, double eddingtonRatio) {
    double parameters = 0.01971853;
    double xi = parameters * eddingtonRatio;

    return (calcSpinupNT(spin) + sHD_min(spin) * xi) / (1 + xi);
}

double calcMagnetization(double spin, double eddingtonRatio) {
    // Fitting parameters from RNC23
    const double param_0 = 1.29;
    const double param_1 = 1.88;

    double mdot_Edd = eddingtonRatio * 0.1 / radiativeEfficiencyThin(spin);

    // Clamp inputs
    if(spin > 0.998) spin = 0.998;
    if(spin < -0.998) spin = -0.998;
    if(eddingtonRatio < 0.0) eddingtonRatio = 0.0;
    if(eddingtonRatio > 1000.0) eddingtonRatio = 1000.0;

    if(mdot_Edd > fEdd_ADAF_threshold) {
        // Thin disk regime: magnetization increases with accretion rate
        double evolutionFactor = pow(eddingtonRatio / param_1, param_0);
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

double etaEM(double spin, double eddingtonRatio) {

    double efficiency = jetEfficiencyFunction(calcMagnetization(spin, eddingtonRatio), spin);
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


double calcJetPower(double mass, double spin, double eddingtonRatio, 
                   double madness, int mad_superEddington, int mad_subEddington) {

    double mdot_Edd = eddingtonRatio * 0.1 / radiativeEfficiencyThin(spin);
    double massInflowRate = calcEddingtonAccretionRateFiducial(mass) * mdot_Edd;
    int jets_enabled = (mad_superEddington && (mdot_Edd > fEdd_ADAF_threshold)) || 
                       (mad_subEddington && (mdot_Edd <= fEdd_ADAF_threshold));

    if(jets_enabled) {
        // Jet power = η_EM × Ṁ × c²
        double eta = etaEM(spin, eddingtonRatio);
        
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

double calcMassEvolutionRate(double mass, double spin, double eddingtonRatio, int coherent,
                            double madness, int mad_superEddington, int mad_subEddington) {

    double mdot_Edd = eddingtonRatio * 0.1 / radiativeEfficiencyThin(spin);
    double massInflowRate = calcEddingtonAccretionRateFiducial(mass) * mdot_Edd;
    double epsilon = radiativeEfficiency(spin, eddingtonRatio);

    // Calculate jet efficiency
    double eta = 0.0;
    int jets_enabled = (mad_superEddington && (mdot_Edd > fEdd_ADAF_threshold)) || 
                       (mad_subEddington && (mdot_Edd <= fEdd_ADAF_threshold));

    if(jets_enabled) {
        // eta = P_jet / (Ṁ × c²)
        double jetPower_watts = calcJetPower(mass, spin, eddingtonRatio, madness, mad_superEddington, mad_subEddington);
        
        // Convert to erg/s
        double jetPower_erg_per_s = jetPower_watts * 1.0e7;
        
        // Mass-energy inflow rate: Ṁ × c² in erg/s
        double massInflowRate_g_per_s = massInflowRate * SOLAR_MASS / SEC_PER_YEAR;
        double massEnergyRate_erg_per_s = massInflowRate_g_per_s * C * C;
        
        eta = jetPower_erg_per_s / massEnergyRate_erg_per_s;
    }

    // Loss factor: fraction of inflowing mass that stays on BH
    double lossFactor = 1.0 - epsilon - eta;

    if(!coherent) {
        // Incoherent: compute same for flipped orientation, then average
        double epsilon_flip = radiativeEfficiency(-spin, eddingtonRatio);
        double r_Edd_flip = eddingtonRatio * 0.1 / radiativeEfficiencyThin(-spin);
        
        double eta_flip = 0.0;
        if(jets_enabled) {
            double jetPower_watts_flip = calcJetPower(mass, -spin, r_Edd_flip, madness, mad_superEddington, mad_subEddington);
            double jetPower_erg_per_s_flip = jetPower_watts_flip * 1.0e7;
            double massInflowRate_g_per_s = massInflowRate * SOLAR_MASS / SEC_PER_YEAR;
            double massEnergyRate_erg_per_s = massInflowRate_g_per_s * C * C;
            eta_flip = jetPower_erg_per_s_flip / massEnergyRate_erg_per_s;
        }
        
        double lossFactor_flip = 1.0 - epsilon_flip - eta_flip;
        
        // Average the two loss factors
        lossFactor = (lossFactor + lossFactor_flip) / 2.0;
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


double calcSpinup(double spin, double eddingtonRatio, double madness, int mad_superEddington, int mad_subEddington) {

    double mdot_Edd = eddingtonRatio * 0.1 / radiativeEfficiencyThin(spin);
    double s;
    
    if(mad_subEddington && (mdot_Edd < fEdd_ADAF_threshold)) {
        // MAD ADAF regime (hot accretion)
        s = calcSpinupMAD(spin) * madness + calcSpinupNT(spin) * (1.0 - madness);
    }
    else if(mad_superEddington && (mdot_Edd >= fEdd_ADAF_threshold)) {
        // Thin disk regime (super-Eddington relative to ADAF threshold)
        
        // Start with hydrodynamic spinup
        s = sHD(spin, eddingtonRatio);
        
        // Add electromagnetic spinup if spin is non-zero
        if(spin != 0.0) {
            double eta = etaEM(spin, eddingtonRatio);
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

double calcSpinEvolutionRate(double mass, double spin, double eddingtonRatio, int coherent,
                            double madness, int mad_superEddington, int mad_subEddington) {

    double mdot_Edd = eddingtonRatio * 0.1 / radiativeEfficiencyThin(spin);
    double s_coherent = calcSpinup(spin, eddingtonRatio, madness, mad_superEddington, mad_subEddington);
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
        double s_incoherent = calcSpinup(-spin, eddingtonRatio, madness, mad_superEddington, mad_subEddington);
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
                          int *substeps_taken) 
{
    double mass = *mass_ptr;
    double spin = initial_spin;
    double time = 0.0;
    double substep = dt_years * 0.01;  // Initial substep (1% of full timestep)
    
    const int MAX_SUBSTEPS = 500;  // Hard limit as requested
    int substep_count = 0;


    // Integrate until we reach dt_years
    while(time < dt_years) {

        if(mass < 1e-10) {
            if(substeps_taken != NULL) *substeps_taken = substep_count;
            *mass_ptr = 1e-10;  // Floor value
            return spin;
        }

        // Early exit: spin saturated at maximum
        if(fabs(spin) >= 0.998) {
            if(substeps_taken != NULL) *substeps_taken = substep_count;
            *mass_ptr = mass;
            return spin;
        }

        if(substep_count % 10 == 0) {
            double da_dt_now = calcSpinEvolutionRate(mass, spin, initial_mdot_Edd, coherent, 
                madness, mad_superEddington, mad_subEddington);
            double remaining_spin_change = fabs(da_dt_now) * (dt_years - time);
            if(remaining_spin_change < 1e-5) {
                if(substeps_taken != NULL) *substeps_taken = substep_count;
                *mass_ptr = mass;
                return spin;
            }
        }

        // if(substep_count >= MAX_SUBSTEPS) {
        //     if(remaining_spin_change > 1e-2 && fabs(spin - initial_spin) > 1e-3) {
        //         static int warn_count = 0;
        //         warn_count++;
        //         // Only print first 3, then summarize every 1000
        //         if(warn_count <= 3 || warn_count % 1000 == 0) {
        //             fprintf(stderr, "WARNING: MAX_SUBSTEPS=%d hit at %.1f%% completion, "
        //                     "spin=%.4f, remaining_change=%.3e. (occurrence #%d)\n",
        //                     MAX_SUBSTEPS, 100.0*time/dt_years, spin, 
        //                     remaining_spin_change, warn_count);
        //         }
        //     }
        //     break;
        // }

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

                // Recalculate mdot_Edd at this prediction point
                double mdot_edd_fiducial = calcEddingtonAccretionRateFiducial(mass_pred);
                double mdot_actual = (initial_mdot_Edd * 0.1 / radiativeEfficiencyThin(initial_spin)) * mdot_edd_fiducial;
                double mdot_Edd_pred = mdot_actual / mdot_edd_fiducial;

                
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
            stepsize_factor_mass = 0.9 * pow(allowedFractionalMassError / (truncation_error_mass / mass), 0.2);
            stepsize_factor_spin = 0.9 * pow(allowedSpinError / truncation_error_spin, 0.2);

            // double stepsize_factor_spin = 0.9 * pow(allowedSpinError / (truncation_error + 1e-30), 0.2);
            
            // Accept step if within tolerance (matches Python line 368)
            if((stepsize_factor_mass >= 1.0 && stepsize_factor_spin >= 1.0) || 
               substep < dt_years * minimumFractionalTimeResolution) {
                step_accepted = 1;
                mass = proposed_mass;
                spin = proposed_spin;
                time += substep;
                substep_count++;
            }
            
            // Adjust timestep (matches Python line 372)
            // Limit increase to factor of 100
            double stepsize_factor = fmin(stepsize_factor_mass, stepsize_factor_spin);
            substep *= fmin(stepsize_factor, 2.0);
            
            // Ensure minimum resolution
            if(minimumFractionalTimeResolution > 0 && substep < dt_years * minimumFractionalTimeResolution) {
                substep = dt_years * minimumFractionalTimeResolution;
            }
        }
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
    double mdot_Edd_final, eddington_ratio, magnetization, jet_power, jet_efficiency, radiative_eff;
    
    mass_msun = Gal[p].BlackHoleMass * 1e10 / Hubble_h;
    spin = Gal[p].BlackHoleSpin;
    
    // // Skip if BH is too small
    // if(mass < 1e-10 || accreted_mass < 1e-20) {
    //     Gal[p].BlackHoleSpin = 0.0;
    //     Gal[p].BH_Magnetization = 0.0; 
    //     Gal[p].BH_JetPower = 0.0;
    //     Gal[p].Eddrate = 0.0;
    //     Gal[p].BH_JetEfficiency = 0.0;  // Set default
    //     Gal[p].BH_thindiskRadiativeEfficiency = 0.0;
    //     return;
    // }
    
    // safety check
    // Clamp spin BEFORE any calculations
    if(spin > 0.998) spin = 0.998;
    if(spin < -0.998) spin = -0.998;
    if(spin != spin) spin = SpinSeedValue;  // NaN check
    
    // safety check

    //if(mdot_Edd < 1e-10) {
    //    Gal[p].Eddrate = 0.0;
    //    eddington_ratio = 0.0;
    //    Gal[p].BH_Magnetization = calcMagnetization(spin, eddington_ratio);
    //    Gal[p].BH_JetPower = 0.0;
    //    Gal[p].BH_JetEfficiency = 0.0;
    //    Gal[p].BH_thindiskRadiativeEfficiency = radiativeEfficiencyThin(spin);
    //    return;
    //}

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
                                   minimumFractionalTimeResolution, &substeps_taken);
    

    // For incoherent: randomly choose prograde or retrograde (matches Python line 381-382)
    if(!coherent) {
        new_spin = (drand48() < 0.5) ? fabs(new_spin) : -fabs(new_spin);
    }
    
    Gal[p].BlackHoleMass = new_mass_msun * Hubble_h / 1e10;
    Gal[p].BlackHoleSpin = new_spin;
 
    double mdot_edd_final = calcEddingtonAccretionRateFiducial(new_mass_msun);
    mdot_Edd_final = mdot_actual / mdot_edd_final;
    if(mdot_Edd_final > 1000.0) mdot_Edd_final = 1000.0;
    
    // Calculate final properties using NEW mdot_Edd
    eddington_ratio = mdot_Edd_final * radiativeEfficiencyThin(new_spin) / 0.1;  // f_Edd
    magnetization = calcMagnetization(new_spin, eddington_ratio);
    jet_power = calcJetPower(new_mass_msun, new_spin, eddington_ratio, madness, mad_superEddington, mad_subEddington);
    jet_efficiency = etaEM(new_spin, eddington_ratio);
    radiative_eff = radiativeEfficiencyThin(new_spin);

    // Store final values
    Gal[p].Eddrate = mdot_Edd_final;  // Use NEW mdot_Edd
    Gal[p].BH_Magnetization = magnetization;
    Gal[p].BH_JetPower = jet_power;
    Gal[p].BH_JetEfficiency = jet_efficiency;
    Gal[p].BH_thindiskRadiativeEfficiency = radiative_eff;
}

void update_black_hole_spin_coherent(int p, double accreted_mass, double dt) {
    update_black_hole_spin(p, accreted_mass, dt, 1, 1.0, 1, 1);  // coherent, madness=1.0, both regimes enabled
}

void update_black_hole_spin_incoherent(int p, double accreted_mass, double dt) {
    update_black_hole_spin(p, accreted_mass, dt, 0, 1.0, 1, 1);  // incoherent, madness=1.0, both regimes enabled
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
