#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#ifdef MPI
#include <mpi.h>
#endif
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void construct_galaxies(int halonr, int tree)
{
  static int halosdone = 0;
  int prog, fofhalo, ngal;

  HaloAux[halonr].DoneFlag = 1;
  halosdone++;

  prog = Halo[halonr].FirstProgenitor;
  while(prog >= 0)
  {
    if(HaloAux[prog].DoneFlag == 0)
      construct_galaxies(prog, tree);
    prog = Halo[prog].NextProgenitor;
  }

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 0)
  {
    HaloAux[fofhalo].HaloFlag = 1;
    while(fofhalo >= 0)
    {
      prog = Halo[fofhalo].FirstProgenitor;
      while(prog >= 0)
      {
        if(HaloAux[prog].DoneFlag == 0)
          construct_galaxies(prog, tree);
        prog = Halo[prog].NextProgenitor;
      }

      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }
  }

  // At this point, the galaxies for all progenitors of this halo have been
  // properly constructed. Also, the galaxies of the progenitors of all other 
  // halos in the same FOF-group have been constructed as well. We can hence go
  // ahead and construct all galaxies for the subhalos in this FOF halo, and
  // evolve them in time. 

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 1)
  {
    ngal = 0;
    HaloAux[fofhalo].HaloFlag = 2;

    while(fofhalo >= 0)
    {
      ngal = join_galaxies_of_progenitors(fofhalo, ngal);
      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }

    evolve_galaxies(Halo[halonr].FirstHaloInFOFgroup, ngal);
  }

}



int join_galaxies_of_progenitors(int halonr, int ngalstart)
{
  int ngal, prog, i, j, first_occupied, centralgal;
    double previousMvir, previousVvir, previousVmax, SpinMag;//, jHotRatio;
  int step;
    long long HaloLen, lenmax, lenoccmax;

  lenmax = 0;
  lenoccmax = 0;
  first_occupied = Halo[halonr].FirstProgenitor;
  prog = Halo[halonr].FirstProgenitor;

  if(prog >=0)
    if(HaloAux[prog].NGalaxies > 0)
    lenoccmax = -1;

  // Find most massive progenitor that contains an actual galaxy
  // Maybe FirstProgenitor never was FirstHaloInFOFGroup and thus got no galaxy

  while(prog >= 0)
  {
      if(Halo[prog].Len > 0)
          HaloLen = (long long) Halo[prog].Len;
      else
          HaloLen = (long long) (dmax(Halo[prog].Mvir1, dmax(Halo[prog].Mvir2, Halo[prog].Mvir3)) / PartMass);
      
    if(HaloLen > lenmax)
    {
      lenmax = HaloLen;
    }
    if(lenoccmax != -1 && HaloLen > lenoccmax && HaloAux[prog].NGalaxies > 0)
    {
      lenoccmax = HaloLen;
      first_occupied = prog;
    }
    prog = Halo[prog].NextProgenitor;
  }

  ngal = ngalstart;
  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0)
  {
    for(i = 0; i < HaloAux[prog].NGalaxies; i++)
    {
      if(ngal >= FoF_MaxGals)
      {
        printf("Opps. We reached the maximum number FoF_MaxGals=%d of galaxies in FoF... exiting.\n", FoF_MaxGals);
        ABORT(1);
      }

      // This is the cruical line in which the properties of the progenitor galaxies 
      // are copied over (as a whole) to the (temporary) galaxies Gal[xxx] in the current snapshot 
      // After updating their properties and evolving them 
      // they are copied to the end of the list of permanent galaxies HaloGal[xxx] 

      Gal[ngal] = HaloGal[HaloAux[prog].FirstGalaxy + i];
      Gal[ngal].HaloNr = halonr;
      Gal[ngal].dT = -1.0;

      // this deals with the central galaxies of (sub)halos 
      if(Gal[ngal].Type == 0 || Gal[ngal].Type == 1)
      {

        // remember properties from the last snapshot
        previousMvir = Gal[ngal].Mvir;
        previousVvir = Gal[ngal].Vvir;
        previousVmax = Gal[ngal].Vmax;
        Gal[ngal].prevRvir = Gal[ngal].Rvir;
          Gal[ngal].prevRhot = sqrt(Gal[ngal].R2_hot_av);
          Gal[ngal].prevVhot = (2 * Gal[ngal].Vvir * Gal[ngal].CoolScaleRadius) / Gal[ngal].prevRhot ;

        if(prog == first_occupied)
        {
          // update properties of this galaxy with physical properties of halo 
          Gal[ngal].MostBoundID = Halo[halonr].MostBoundID;

          for(j = 0; j < 3; j++)
          {
            Gal[ngal].Pos[j] = Halo[halonr].Pos[j];
            Gal[ngal].Vel[j] = Halo[halonr].Vel[j];
          }
  
            
            if(Halo[halonr].Len > 0)
                HaloLen = (long long) Halo[halonr].Len;
            else
                HaloLen = (long long) (dmax(Halo[halonr].Mvir1, dmax(Halo[halonr].Mvir2, Halo[halonr].Mvir3)) / PartMass);
            
          if(HaloLen > Gal[ngal].LenMax) Gal[ngal].LenMax = HaloLen;
          Gal[ngal].Len = HaloLen;
          Gal[ngal].Vmax = Halo[halonr].Vmax;
          assert(Gal[ngal].LenMax>=Gal[ngal].Len);
            
              
            if(halonr == Halo[halonr].FirstHaloInFOFgroup)
            {
                Gal[ngal].deltaMvir = get_virial_mass(halonr, ngal) - Gal[ngal].Mvir;
                Gal[ngal].Vvir = get_virial_velocity(halonr, ngal); // Keeping Vvir fixed for subhaloes
                Gal[ngal].Mratio = Gal[ngal].Mvir / (Gal[ngal].Len * PartMass);
                Gal[ngal].Mvir = get_virial_mass(halonr, ngal);
                Gal[ngal].Rvir = get_virial_radius(halonr, ngal, Gal[ngal].Mvir);
            }
            else
            {
                double Mvir_sat = get_virial_mass(halonr, ngal);
                if(Gal[ngal].Mratio > 0.0 && Gal[ngal].Mratio < 1.0) Mvir_sat *= Gal[ngal].Mratio;
                Gal[ngal].deltaMvir = Mvir_sat - Gal[ngal].Mvir;
                Gal[ngal].Mvir = Mvir_sat;
                Gal[ngal].Rvir = get_virial_radius(halonr, ngal, Mvir_sat);
            }
                
        Gal[ngal].DiskScaleRadius = get_disk_radius(halonr, ngal); // Only update the scale radius for centrals.  Satellites' spin will be too variable and untrustworthy for new cooling.

          Gal[ngal].Cooling = 0.0;
          Gal[ngal].Heating = 0.0;
          Gal[ngal].SNreheatRate = 0.0;
          Gal[ngal].SNejectRate = 0.0;

          for(step = 0; step < STEPS; step++)
          {
//            Gal[ngal].SfrDisk[step] = Gal[ngal].SfrBulge[step]  = 0.0;
            Gal[ngal].SfrFromH2[step] = Gal[ngal].SfrInstab[step] = Gal[ngal].SfrMerge[step] = 0.0;
            Gal[ngal].SfrDiskColdGas[step] = Gal[ngal].SfrDiskColdGasMetals[step] = 0.0;
            Gal[ngal].SfrBulgeColdGas[step] = Gal[ngal].SfrBulgeColdGasMetals[step] = 0.0;
          }

            
          SpinMag = sqrt(sqr(Halo[halonr].Spin[0]) + sqr(Halo[halonr].Spin[1]) + sqr(Halo[halonr].Spin[2]));
          if(SpinMag>0)
          {
            for(j = 0; j < 3; j++) // Also only update the spin direction of the hot gas for centrals
                Gal[ngal].SpinHot[j] = Halo[halonr].Spin[j] / SpinMag;
          }

          if(halonr == Halo[halonr].FirstHaloInFOFgroup)
          {
            // a central galaxy
            Gal[ngal].mergeType = 0;
            Gal[ngal].mergeIntoID = -1;
            Gal[ngal].MergTime = 999.9;            

            Gal[ngal].Type = 0;
          }
          else
          {
            // a satellite with subhalo
              Gal[ngal].mergeType = 0;
              Gal[ngal].mergeIntoID = -1;
              
            if(Gal[ngal].Type == 0)  // remember the infall properties before becoming a subhalo
            {
              Gal[ngal].infallMvir = previousMvir;
              Gal[ngal].infallVvir = previousVvir;
              Gal[ngal].infallVmax = previousVmax;
            }

            Gal[ngal].Type = 1;
          }
        }
        else
        { // NO ORPHANS SHOULD BE PRESENT IN THE MODEL
          // an orhpan satellite galaxy
          Gal[ngal].deltaMvir = -1.0*Gal[ngal].Mvir;
          Gal[ngal].Mvir = 0.0;

          if(Gal[ngal].MergTime > 999.0 || Gal[ngal].Type == 0)
          {
            // Here the galaxy has gone from type 0 to type 2. Merge it!
            Gal[ngal].MergTime = 0.0;
          
            Gal[ngal].infallMvir = previousMvir;
            Gal[ngal].infallVvir = previousVvir;
            Gal[ngal].infallVmax = previousVmax;
          }

          Gal[ngal].Type = 2;
        }
      }
        
        if(Gal[ngal].TypeMax < Gal[ngal].Type) Gal[ngal].TypeMax = Gal[ngal].Type;

      // Note: Galaxies that are already type 2 do not need special treatment at this point 
      if(Gal[ngal].Type < 0 || Gal[ngal].Type > 2)
      {
        printf("what's that????\n");
        ABORT(88);
      }

      ngal++;

    }

    prog = Halo[prog].NextProgenitor;
  }

  if(ngal == 0)
  {
    // We have no progenitors with galaxies. This means we create a new galaxy. 
    init_galaxy(ngal, halonr);
    ngal++;
  }

  // Per Halo there can be only one Type 0 or 1 galaxy, all others are Type 2 
  // In fact this galaxy is very likely to be the first galaxy in the halo if first_occupied==FirstProgenitor 
  // and the Type0/1 galaxy in FirstProgenitor was also the first one 
  // This cannot be guaranteed though for the pathological first_occupied != FirstProgenitor case 

  for(i = ngalstart, centralgal = -1; i < ngal; i++)
  {
    if(Gal[i].Type == 0 || Gal[i].Type == 1)
    {
      if(centralgal != -1)
      {
        printf("can't be\n");
        ABORT(8);
      }
      centralgal = i;
    }
  }

  for(i = ngalstart; i < ngal; i++)
    Gal[i].CentralGal = centralgal;

  return ngal;

}



void evolve_galaxies(int halonr, int ngal)	// note: halonr is here the FOF-background subhalo (i.e. main halo)
{
  int p, i, step, centralgal, merger_centralgal, currenthalo, offset, k_now;
  double infallingGas, coolingGas, deltaT, time, galaxyBaryons, currentMvir, DiscGasSum, dt, InstantTimeFrame, StellarOutput[2];
  int gas_instab = 0;

  centralgal = Gal[0].CentralGal;
  if(Gal[centralgal].Type != 0 || Gal[centralgal].HaloNr != halonr)
  {
    printf("Something wrong here ..... \n");
    ABORT(54);
  }

    assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
  // calculate cosmological gas to be added to halo to maintain baryon fraction
  infallingGas = infall_recipe(centralgal, ngal, ZZ[Halo[halonr].SnapNum]);

  // move mass in outflow and ejected reservoirs based on halo growth and update time-scales accordingly
    assert(Gal[centralgal].OutflowGas >= 0.0);
  if(infallingGas > 0)
      update_galactic_fountain_from_growth(centralgal);
  assert(Gal[centralgal].OutflowGas >= 0.0);
  
  for(p = 0; p < ngal; p++)
  {
      // Reset SFRs for each galaxy
      for(i=0; i<N_BINS; i++)
        Gal[p].DiscSFR[i] = 0.0;

      // check for freshly infallen satellites that need a merger timescale calculated
      if(Gal[p].Type == 1 && Gal[p].MergTime > 999.0)
        Gal[p].MergTime = estimate_merging_time(p, centralgal);
  }

  // we integrate things forward by using a number of intervals equal to STEPS 
  for(step = 0; step < STEPS; step++)
  {

    // Loop over all galaxies in the halo 
    for(p = 0; p < ngal; p++)
    {
	  DiscGasSum = get_disc_gas(p);
	  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
	  assert(Gal[p].MetalsColdGas <= Gal[p].ColdGas);
          assert(Gal[p].MetalsColdGas <= DiscGasSum);

	  if(step==0 || HotStripOn<3) Gal[p].MaxStrippedGas = Gal[p].HotGas / STEPS;
        
      // don't treat galaxies that have already merged 
      if(Gal[p].mergeType > 0)
        continue;

      // Debugging
      assert(Gal[p].MetalsHotGas>=0);


      deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
      dt = (deltaT / STEPS);
      time = Age[Gal[p].SnapNum] - (step + 0.5) * dt;
      k_now = get_stellar_age_bin_index(time);

      if(Gal[p].dT < 0.0)
        Gal[p].dT = deltaT;
    
      // hot-gas stripping of satellites --  do before reincorporation to ensure stripping is preventative
        assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
      if(HotStripOn>0 && Gal[p].Type == 1 && Gal[p].HotGas > 0.0 && Gal[p].MaxStrippedGas>0.0)

	    // Debugging
	    assert(Gal[p].MetalsHotGas>=0);
	    // Debugging

            Gal[p].MaxStrippedGas = strip_from_satellite(halonr, centralgal, p, Gal[p].MaxStrippedGas, k_now);
            assert(Gal[p].MetalsHotGas>=0);
            assert(Gal[p].MetalsHotGas<=Gal[p].HotGas);
            assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);

      // fresh cosmological added to CGM
      if(p==centralgal)
      {
          assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
          add_infall_to_hot(p, infallingGas / STEPS);
          assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
          assert(Gal[p].MetalsHotGas>=0);
          assert(Gal[p].MetalsHotGas<=Gal[p].HotGas);
      }
        
      assert(Gal[p].OutflowGas >= 0.0);

      // move outflowing, fountaining, and reincorporating ejected gas as appropriate
      reincorporate_gas(p, centralgal, dt);
      assert(Gal[p].MetalsHotGas<=Gal[p].HotGas);
      assert(Gal[p].OutflowGas >= 0.0);
          
      // Ram pressure stripping of cold gas from satellites
      if(RamPressureOn>0 && Gal[p].Type == 1 && Gal[p].ColdGas>0.0)
          ram_pressure_stripping(centralgal, p, k_now);
      assert(Gal[p].MetalsHotGas<=Gal[p].HotGas);
      assert(Gal[p].OutflowGas >= 0.0);

      // -----------------------------------------------------------------------
      // PHASE 1: compute raw cooling gas and accumulate BH supply reservoirs
      // -----------------------------------------------------------------------
      //
      // cooling_recipe_part_a() replaces the old cooling_recipe() + do_AGN_heating() pair.
      // It returns the raw cooling gas *before* any AGN suppression and outputs the
      // specific energy change needed by bh_feedback() later.  It also computes the
      // radio-mode AGN accretion rate and deposits it into BH_incoherent_reservoir,
      // deducting the corresponding mass (and metals) from HotGas.
      //
      // check_disk_instability() deposits any cold-disc BH supply into BH_coherent_reservoir
      // (via deal_with_unstable_gas) rather than evolving the spin directly.  It must run
      // before cool_gas_onto_galaxy() so the disc state it sees is self-consistent (stale-disc
      // approximation: the error shrinks with more sub-steps, STEPS).
      //
      // Preventing cooling when there's no angular momentum, as the code isn't built to
      // handle that.  This only cropped up for Vishnu haloes with 2 particles, which clearly
      // weren't interesting/physical.  Haloes always have some spin otherwise.
      double coolingGas_raw      = 0.0;
      double specific_energy_change = 0.0;
      int    do_cooling = !(Gal[p].SpinHot[0]==0 && Gal[p].SpinHot[1]==0 && Gal[p].SpinHot[2]==0);

      if(do_cooling)
          coolingGas_raw = cooling_recipe_part_a(p, dt, &specific_energy_change);

      // Check for disk instability (before cool_gas_onto_galaxy to use stale disc state).
      // deal_with_unstable_gas() adds cold-disc BH supply to BH_coherent_reservoir.
      if(DiskInstabilityOn>0)
          gas_instab = check_disk_instability(p, dt, step, time, k_now);

      // -----------------------------------------------------------------------
      // PHASE 2: single unified BH evolution
      // -----------------------------------------------------------------------
      //
      // Capture reservoir values before evolve_black_hole() zeroes them; bh_feedback()
      // needs the original supply amounts to compute heating and quasar winds.
      double radio_supply  = Gal[p].BH_incoherent_reservoir;
      double instab_supply = Gal[p].BH_coherent_reservoir;

      // One spin update per sub-step using both reservoirs.  Returns the radio-mode
      // net BH mass increment (after losses), needed by bh_feedback() for jet heating.
      double radio_net = evolve_black_hole(p, dt);

      // DeltaBlackHoleMass now holds the instability-driven fraction of the net increment.
      // Accumulate it into InstaBlackHoleMass and record the per-sub-step rate.
      Gal[p].InstaBlackHoleMass += Gal[p].DeltaBlackHoleMass;
      Gal[p].InstaBHaccreted     = Gal[p].DeltaBlackHoleMass / dt;

      // -----------------------------------------------------------------------
      // PHASE 3: AGN feedback and gas cooling onto the disc
      // -----------------------------------------------------------------------
      //
      // bh_feedback() computes radio-mode jet heating (suppressing coolingGas) and
      // fires the quasar-mode wind.  It must run after evolve_black_hole() so that
      // BH_JetPower (for SpinEvolutionOn>0) and BH_thindiskRadiativeEfficiency (for
      // quasar_mode_wind) both reflect the combined spin update just performed.
      coolingGas = bh_feedback(p, radio_supply, instab_supply, radio_net,
                               coolingGas_raw, specific_energy_change, dt);

      if(do_cooling)
      {
          // this is no longer an accurate representation of the actual energy change in the
          // gas as it cools, but we keep it for consistency with the old accounting.
          if(coolingGas > 0.0)
              Gal[p].Cooling += coolingGas * specific_energy_change;

          Gal[p].AccretedGasMass += coolingGas;
          cool_gas_onto_galaxy(p, coolingGas);
      }

      // debugging
      if(p == centralgal) assert(Gal[centralgal].MetalsLocalIGS >= 0);
      // debugging


      // Update radii of the annuli.  Done again at end of full time-step, so seems redundant on first sub-time-step (costly)
      assert(Gal[p].OutflowGas >= 0.0);
      if(Gal[p].Mvir > 0 && Gal[p].Rvir > 0 && step!=0) update_disc_radii(p);
      assert(Gal[p].OutflowGas >= 0.0);

      // If using newer feedback model, calculate the instantaneous recycling fraction for the current time-step + apply delayed feedback from earlier stellar populations
      if(DelayedFeedbackOn>0 && N_AGE_BINS>1 && (Gal[p].StellarMass>0 || Gal[p].ICS>0 || Gal[p].LocalIGS>0))
      {
        InstantTimeFrame = 0.5*(AgeBinEdge[k_now+1] - AgeBinEdge[k_now]);
        get_RecycleFraction_and_NumSNperMass(0.0, InstantTimeFrame, StellarOutput);
        RecycleFraction = 1.0*StellarOutput[0];
        SNperMassFormed = 1.0*StellarOutput[1];

        assert(Gal[p].OutflowGas >= 0.0);
        delayed_feedback(p, k_now, time, dt);
        assert(Gal[p].OutflowGas >= 0.0);
      }

      // debugging
      if(p == centralgal) assert(Gal[centralgal].MetalsLocalIGS >= 0);
      // debugging

      // precess gas disc
      if(GasPrecessionOn && Gal[p].StellarMass>0.0 && get_disc_gas(p)>0.0)
        precess_gas(p, dt);
      assert(Gal[p].SpinGas[0]==Gal[p].SpinGas[0]);

      // stars form and then explode!
      starformation_and_feedback(p, centralgal, dt, step, time, k_now);
      assert(Gal[p].OutflowGas >= 0.0);

      // debugging
      if(p == centralgal) assert(Gal[centralgal].MetalsLocalIGS >= 0);
      // debugging

    }

    // check for satellite disruption and merger events 
    for(p = 0; p < ngal; p++)
    {

      if((Gal[p].Type == 1 || Gal[p].Type == 2) && Gal[p].mergeType == 0)  // satellite galaxy!
      {
        if(Gal[p].MergTime > 999.0)
        {
          printf("satellite doesn't have a merging time! %f\n", Gal[p].MergTime);
          //ABORT(77);
        }

        deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
        dt = (deltaT / STEPS);
        Gal[p].MergTime -= dt;
        
        // only consider mergers or disruption for halo-to-baryonic mass ratios below the threshold
        // or for satellites with no baryonic mass (they don't grow and will otherwise hang around forever)
        currentMvir = Gal[p].Mvir - Gal[p].deltaMvir * (1.0 - ((double)step + 1.0) / (double)STEPS);
        galaxyBaryons = Gal[p].StellarMass + Gal[p].ColdGas + Gal[p].HotGas + Gal[p].BlackHoleMass + Gal[p].FountainGas + Gal[p].OutflowGas + Gal[p].ICS + Gal[p].ICBHmass; // not adding ejected gas as it outside the virial radius, strictly
          
        double reionization_modifier;
        if(ReionizationOn)
            reionization_modifier = do_reionization(p, ZZ[Halo[halonr].SnapNum]);
        else
            reionization_modifier = 1.0;
          
        if((galaxyBaryons == 0.0) || (galaxyBaryons > 0.0 && (currentMvir / galaxyBaryons <= ThresholdSatDisruption / reionization_modifier)))        
        {
          if(Gal[p].Type==1) 
            merger_centralgal = centralgal;
          else
            merger_centralgal = Gal[p].CentralGal;

          if(Gal[merger_centralgal].mergeType > 0) 
            merger_centralgal = Gal[merger_centralgal].CentralGal;

          Gal[p].mergeIntoID = NumGals + merger_centralgal;  // position in output 
            
          time = Age[Gal[p].SnapNum] - (step + 0.5) * dt;
          k_now = get_stellar_age_bin_index(time);
            
          if(Gal[p].MergTime > 0.0)  // disruption has occured!
          {
              assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
            disrupt_satellite_to_ICS(centralgal, p, k_now);
              assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
          }
          else
          {
              assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
            deal_with_galaxy_merger(p, merger_centralgal, centralgal, time, dt, step, k_now);
              assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);
          }
 
        }
     
      }

      // debuggin lines
      if(Gal[centralgal].MetalsLocalIGS < 0) {
          printf("ERROR after processing galaxy %d (Type=%d)\n", p, Gal[p].Type);
          printf("  MetalsLocalIGS = %.6e\n", Gal[centralgal].MetalsLocalIGS);
      }
      assert(Gal[centralgal].MetalsLocalIGS >= 0);
      // debuggin lines

    }

    // debuggin lines
    if(Gal[centralgal].MetalsLocalIGS < 0 || Gal[centralgal].MetalsLocalIGS != Gal[centralgal].MetalsLocalIGS) {
        printf("\n=== METAL CONSERVATION ERROR ===\n");
        printf("STEP: %d/%d\n", step, STEPS);
        printf("centralgal: %d\n", centralgal);
        printf("MetalsLocalIGS: %.10e\n", Gal[centralgal].MetalsLocalIGS);
        printf("LocalIGS: %.10e\n", Gal[centralgal].LocalIGS);
        printf("\nSatellites in this halo:\n");

        for(int i = 0; i < ngal; i++) {
            if(i == centralgal) continue;
            if(Gal[i].Type == 1 || Gal[i].Type == 2) {
                printf("  Galaxy %d: Type=%d, mergeType=%d, MergTime=%.2f\n", i, Gal[i].Type, Gal[i].mergeType, Gal[i].MergTime);
                printf("    LocalIGS=%.3e, MetalsLocalIGS=%.3e\n", Gal[i].LocalIGS, Gal[i].MetalsLocalIGS);
            }
        }
        printf("================================\n");
        fflush(stdout);
      }
     // debuggin lines
      assert(Gal[centralgal].MetalsLocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS != INFINITY);

  } // end move forward in interval STEPS 


  // extra miscellaneous stuff 
  for(p = 0; p < ngal; p++)
  {
    deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
    Gal[p].Cooling /= deltaT;
    Gal[p].Heating /= deltaT;

    Gal[p].BHaccreted = (Gal[p].RadioBHaccreted +  Gal[p].InstaBHaccreted);

    Gal[p].SNreheatRate /= deltaT;
    Gal[p].SNejectRate /= deltaT;

    double solarmass_per_year = 6.31e25;
    Gal[p].Lbol = radiativeEfficiencyThin(Gal[p].BlackHoleSpin)*Gal[p].BHaccreted*solarmass_per_year*(pow(C,2));
      
    assert(Gal[p].SNreheatRate >= 0);
    assert(Gal[p].SNejectRate >= 0);
      
    if(Gal[p].Mvir > 0 && Gal[p].Rvir > 0)
    {
        update_disc_radii(p);
        update_HI_H2(p, time, k_now);
    }
      
    Gal[p].prevHotGasPotential = Gal[p].HotGasPotential;
      
      
  }

  // attach final galaxy list to halo 
  offset = 0;
  for(p = 0, currenthalo = -1; p < ngal; p++)
  {
    
    if(Gal[p].HaloNr != currenthalo)
    {
      currenthalo = Gal[p].HaloNr;
      HaloAux[currenthalo].FirstGalaxy = NumGals;
      assert(NumGals>=0);
      HaloAux[currenthalo].NGalaxies = 0;
    }
      
    // Merged galaxies won't be output. So go back through its history and find it
    // in the previous timestep. Then copy the current merger info there.
    offset = 0;
    i = p-1;
    while(i >= 0)
    {
     if(Gal[i].mergeType > 0) 
       if(Gal[p].mergeIntoID > Gal[i].mergeIntoID)
         offset++;  // these galaxies won't be kept so offset mergeIntoID below
     i--;
    }
    
    i = -1;
    if(Gal[p].mergeType > 0)
    {
      i = NumGals;
      while(i >= 0)
      {
        if(HaloGal[i].GalaxyNr == Gal[p].GalaxyNr)
          break;
        else
          i--;
      }
      
      if(i < 0)
      {
        printf("Ran over the end of HaloGal looking for progenitor!\n");
        ABORT(12);
      }
      
      HaloGal[i].mergeType = Gal[p].mergeType;
      HaloGal[i].mergeIntoID = Gal[p].mergeIntoID - offset;
      HaloGal[i].mergeIntoGalaxyNr = Gal[p].mergeIntoGalaxyNr;
      HaloGal[i].mergeIntoSnapNum = Halo[currenthalo].SnapNum;
    }
    
    if(Gal[p].mergeType == 0)
    {
      if(NumGals >= MaxGals)
      {
        printf("maximum number of galaxies reached... exiting.\n");
        ABORT(1);
      }

      Gal[p].SnapNum = Halo[currenthalo].SnapNum;
      HaloGal[NumGals++] = Gal[p];
      HaloAux[currenthalo].NGalaxies++;
    }
  }
}






