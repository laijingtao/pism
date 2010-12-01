// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include "iceModel.hh"

#include <petscsys.h>
#include <stdarg.h>
#include <stdlib.h>

//!  Computes volume and area of ice sheet, for reporting purposes.
/*!
Also computes ice volumes in the three mask regions SHEET, DRAGGING, FLOATING.

Communication done for global max and global sum.

Returns area in units of m^2 and volume in m^3.
 */
PetscErrorCode IceModel::volumeArea(PetscScalar& gvolume, PetscScalar& garea,
                                    PetscScalar& gvolSIA, PetscScalar& gvolstream, 
                                    PetscScalar& gvolshelf) {

  PetscErrorCode  ierr;
  PetscScalar     **H;
  PetscScalar     volume=0.0, area=0.0, volSIA=0.0, volstream=0.0, volshelf=0.0;
  
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  // FIXME: we should use cell_area instead
  const PetscScalar   a = grid.dx * grid.dy; // area unit (m^2)
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        area += a;
        const PetscScalar dv = a * H[i][j];
        volume += dv;
        if (vMask.value(i,j) == MASK_SHEET)   volSIA += dv;
        else if (vMask.value(i,j) == MASK_DRAGGING_SHEET)   volstream += dv;
        else if (vMask.is_floating(i,j))   volshelf += dv;
      }
    }
  }  
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volume, &gvolume, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&area, &garea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volSIA, &gvolSIA, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volstream, &gvolstream, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&volshelf, &gvolshelf, grid.com); CHKERRQ(ierr);
  return 0;
}


/*!
Computes fraction of the base which is melted and the ice basal temperature at the
center of the ice sheet.

Communication occurs here.
 */
PetscErrorCode IceModel::energyStats(PetscScalar iarea, PetscScalar &gmeltfrac,
                                     PetscScalar &gtemp0) {
  PetscErrorCode    ierr;
  PetscScalar       meltarea = 0.0, temp0 = 0.0;
  const PetscScalar a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  
  ierr = vH.begin_access(); CHKERRQ(ierr);

  // use Enth3 to get stats
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
	// accumulate area of base which is at melt point
	if (EC->isTemperate(Enthbase[i][j], EC->getPressureFromDepth(vH(i,j)) ))  
	  meltarea += a;
      }
      // if you happen to be at center, record absolute basal temp there
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
	ierr = EC->getAbsTemp(Enthbase[i][j],EC->getPressureFromDepth(vH(i,j)), temp0);
	CHKERRQ(ierr);
      }
    }
  }
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  ierr = vH.end_access(); CHKERRQ(ierr);

  // communication
  ierr = PetscGlobalSum(&meltarea, &gmeltfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&temp0,    &gtemp0,    grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
  if (iarea > 0.0)   gmeltfrac = gmeltfrac / iarea;
  else gmeltfrac = 0.0;

  return 0;
}


/*!
Computes fraction of the ice which is as old as the start of the run (original).
Communication occurs here.
 */
PetscErrorCode IceModel::ageStats(PetscScalar ivol, PetscScalar &gorigfrac) {
  PetscErrorCode  ierr;

  gorigfrac = -1.0;  // result value if not do_age

  if (!config.get_flag("do_age")) 
    return 0;  // leave now

  const PetscScalar  a = grid.dx * grid.dy * 1e-3 * 1e-3, // area unit (km^2)
                     currtime = grid.year * secpera; // seconds

  PetscScalar **H, *tau, origvol = 0.0;
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);

  // compute local original volume
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        // accumulate volume of ice which is original
        ierr = tau3.getInternalColumn(i,j,&tau); CHKERRQ(ierr);
        const PetscInt  ks = grid.kBelowHeight(H[i][j]);
        for (PetscInt k=1; k<=ks; k++) {
          // ice in segment is original if it is as old as one year less than current time
          if (0.5*(tau[k-1]+tau[k]) > currtime - secpera)
            origvol += a * 1.0e-3 * (grid.zlevels[k] - grid.zlevels[k-1]);
        }
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access(); CHKERRQ(ierr);

  // communicate to turn into global original fraction
  ierr = PetscGlobalSum(&origvol,  &gorigfrac, grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
  if (ivol > 0.0)    gorigfrac = gorigfrac / ivol;
  else gorigfrac = 0.0;

  return 0;
}


PetscErrorCode IceModel::summary(bool tempAndAge) {
  PetscErrorCode  ierr;
  PetscScalar     **H;
  PetscScalar     divideH;
  PetscScalar     gdivideH, gdivideT, gvolume, garea;
  PetscScalar     gvolSIA, gvolstream, gvolshelf;
  PetscScalar     meltfrac = 0.0, origfrac = 0.0;

  // get volumes in m^3 and areas in m^2
  ierr = volumeArea(gvolume, garea, gvolSIA, gvolstream, gvolshelf); CHKERRQ(ierr);
  
  // get thick0 = gdivideH
  ierr = vH.get_array(H); CHKERRQ(ierr);
  divideH = 0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        divideH = H[i][j];
      }
    }
  }  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&divideH, &gdivideH, grid.com); CHKERRQ(ierr);

  if (tempAndAge || (getVerbosityLevel() >= 3)) {
    ierr = energyStats(garea, meltfrac, gdivideT); CHKERRQ(ierr);
  }

  if ((tempAndAge || (getVerbosityLevel() >= 3)) && (config.get_flag("do_age"))) {
    ierr = ageStats(gvolume, origfrac); CHKERRQ(ierr);
  }

  // report CFL violations
  if (CFLviolcount > 0.0) {
    const PetscScalar CFLviolpercent = 100.0 * CFLviolcount / (grid.Mx * grid.Mz * grid.Mz);
    // at default verbosity level, only report CFL viols if above:
    const PetscScalar CFLVIOL_REPORT_VERB2_PERCENT = 0.1;
    if (   (CFLviolpercent > CFLVIOL_REPORT_VERB2_PERCENT) 
        || (getVerbosityLevel() > 2) ) {
      char tempstr[90] = "";
      snprintf(tempstr,90,
              "  [!CFL#=%1.0f (=%5.2f%% of 3D grid)] ",
              CFLviolcount,CFLviolpercent);
      stdout_flags = tempstr + stdout_flags;
    }
  }
   
  // main report: 'S' line
  ierr = summaryPrintLine(PETSC_FALSE,(PetscTruth)tempAndAge,grid.year,dt,
                          gvolume,garea,meltfrac,gdivideH,gdivideT); CHKERRQ(ierr);

  return 0;
}


//! Print a line to stdout which summarizes the state of the modeled ice sheet at the end of the time step.
/*!
Generally, two lines are printed to stdout, the first starting with a space 
and the second starting with the character 'S' in the left-most column (column 1).

The first line shows flags for which processes executed, and the length of the
time step (and/or substeps under option -skip).  See IceModel::run()
for meaning of these flags.

If IceModel::printPrototype is TRUE then the first line does not appear and
the second line has alternate appearance.  Specifically, different column 1
characters are printed:
  - 'P' line gives names of the quantities reported in the 'S' line, the
    "prototype", while
  - 'U' line gives units of these quantities.
This column 1 convention allows automatic tools to read PISM stdout
and produce time-series.  The 'P' and 'U' lines are intended to appear once at
the beginning of the run, while an 'S' line appears at every time step.
This base class version gives a report based on the information included in the 
EISMINT II intercomparison of ice sheet models[\ref EISMINT00].

Note that the inputs \c volume and \c area to this method are in m^3 and m^2,
respectively.  Thus all inputs to this method are in MKS except for \c year.

The resulting numbers on an 'S' line have the following meaning in this base
class version:
  - \c ivol is the total ice sheet volume
  - \c iarea is the total area occupied by positive thickness ice
  - \c thick0 is the ice thickness at the center of the computational domain
  - \c temp0 is the ice basal temperature at the center of the computational domain
The last two can be interpreted as "sanity checks", because they give
information about a location which may or may not be "typical".

For more description and examples, see the PISM User's Manual.
Derived classes of IceModel may redefine this method and print alternate
information.  Use of DiagnosticTimeseries may be superior, however.
 */
PetscErrorCode IceModel::summaryPrintLine(
     PetscTruth printPrototype,  bool tempAndAge,
     PetscScalar year,  PetscScalar delta_t,
     PetscScalar volume,  PetscScalar area,
     PetscScalar /* meltfrac */,  PetscScalar H0,  PetscScalar T0) {

  PetscErrorCode ierr;
  const bool do_temp = config.get_flag("do_temp");

  const int log10scale = static_cast<int>(config.get("summary_volarea_scale_factor_log10"));
  const double scale = pow(10.0, static_cast<double>(log10scale));
  char  volscalestr[10] = "     ", // for special case: blank when 10^0 = 1 scaling
        areascalestr[10] = "   ";  // ditto
  if (log10scale != 0) {
    snprintf(volscalestr, sizeof(volscalestr), "10^%1d_", log10scale);
    strcpy(areascalestr,volscalestr);
  }

  // this version keeps track of what has been done so as to minimize stdout:
  static string stdout_flags_count0;
  static int    mass_cont_sub_counter = 0;  
  static double mass_cont_sub_dtsum = 0.0;
  
  if (printPrototype == PETSC_TRUE) {
    if (do_temp) {
      ierr = verbPrintf(2,grid.com,
          "P         YEAR:     ivol   iarea     thick0     temp0\n");
      ierr = verbPrintf(2,grid.com,
          "U        years %skm^3 %skm^2        m         K\n",
          volscalestr,areascalestr);
    } else {
      ierr = verbPrintf(2,grid.com,
          "P         YEAR:     ivol   iarea     thick0\n");
      ierr = verbPrintf(2,grid.com,
          "U        years %skm^3 %skm^2        m\n",
          volscalestr,areascalestr);
    }
  } else {
    if (mass_cont_sub_counter == 0)
      stdout_flags_count0 = stdout_flags;
    if (delta_t > 0.0) {
      mass_cont_sub_counter++;      
      mass_cont_sub_dtsum += delta_t;
    }
    if ((tempAndAge == PETSC_TRUE) || (!do_temp) || (getVerbosityLevel() > 2)) {
      char tempstr[90] = "";
      const PetscScalar major_dt_years = mass_cont_sub_dtsum / secpera;
      if (mass_cont_sub_counter == 1) {
        snprintf(tempstr,90, " (dt=%.5f)", major_dt_years);
      } else {
        snprintf(tempstr,90, " (dt=%.5f in %d substeps; av dt_sub_mass_cont=%.5f)",
          major_dt_years, mass_cont_sub_counter, major_dt_years / mass_cont_sub_counter);
      }
      stdout_flags_count0 += tempstr;
      if (delta_t > 0.0) { // avoids printing an empty line if we have not done anything
        stdout_flags_count0 += "\n";
        ierr = verbPrintf(2,grid.com, stdout_flags_count0.c_str()); CHKERRQ(ierr);
      }
      if (stdout_ssa.length() > 0) {
        stdout_ssa += "\n";
        ierr = verbPrintf(2,grid.com, stdout_ssa.c_str()); CHKERRQ(ierr);
      }
      if (do_temp) {
        ierr = verbPrintf(2,grid.com, 
          "S %12.5f: %8.5f %7.4f %10.3f %9.4f\n",
          year, volume/(scale*1.0e9), area/(scale*1.0e6), H0, T0); CHKERRQ(ierr);
      } else {
        ierr = verbPrintf(2,grid.com, 
          "S %12.5f: %8.5f %7.4f %10.3f\n",
          year, volume/(scale*1.0e9), area/(scale*1.0e6), H0); CHKERRQ(ierr);
      }
      mass_cont_sub_counter = 0;      
      mass_cont_sub_dtsum = 0.0;
    }
  }
  return 0;
}

//! \brief Sets entrues of result to corresponding processor ranks.
PetscErrorCode IceModel::compute_proc_ice_area(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscInt ice_filled_cells = 0;

  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      if (vH(i,j) > 0) {
	ice_filled_cells += 1;
      }
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      result(i,j) = ice_filled_cells;
  ierr = result.end_access();

  ierr = result.set_name("proc_ice_area"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic",
			  "number of cells containing ice in a processor's domain",
			  "", ""); CHKERRQ(ierr);
  result.time_independent = true;
  return 0;
}


//! \brief Computes f(|v|) as described in [\ref BBssasliding] (page 7, equation 22). 
static PetscScalar bueler_brown_f(PetscScalar v_squared) {
  const PetscScalar inC_fofv = 1.0e-4 * PetscSqr(secpera),
    outC_fofv = 2.0 / pi;
  
  return 1.0 - outC_fofv * atan(inC_fofv * v_squared);
}

//! \brief Computes the map of f(|v|) (see [\ref BBssasliding], equation 22)
PetscErrorCode IceModel::compute_bueler_brown_f(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal fill_value = -0.01;
  IceModelVec2V *vel_ssa;

  ierr = stress_balance->get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);

  ierr = vel_ssa->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = bueler_brown_f((*vel_ssa)(i,j).magnitude_squared());
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vel_ssa->end_access(); CHKERRQ(ierr);

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);

  ierr = result.set_name("bueler_brown_f"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "f(|v|) in Bueler and Brown (2009), equation 22", "", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


//! Compute the CTS field, CTS = E/E_s(p) and put in a global IceModelVec3 provided by user.
PetscErrorCode IceModel::compute_cts(IceModelVec3 &useForCTS) {
  PetscErrorCode ierr;

  ierr = setCTSFromEnthalpy(useForCTS); CHKERRQ(ierr);
  ierr = useForCTS.set_name("cts"); CHKERRQ(ierr);
  ierr = useForCTS.set_attrs(
     "diagnostic",
     "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
     "", ""); CHKERRQ(ierr);
  return 0;
}


//! Compute the liquid fraction, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceModel::compute_liqfrac(IceModelVec3 &useForLiqfrac) {
  PetscErrorCode ierr;

  if (config.get_flag("do_cold_ice_methods")) {
    ierr = useForLiqfrac.set(0.0); CHKERRQ(ierr);
  } else {
    ierr = setLiquidFracFromEnthalpy(useForLiqfrac); CHKERRQ(ierr);
  }
  ierr = useForLiqfrac.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = useForLiqfrac.set_attrs(
     "diagnostic","liquid water fraction in ice (between 0 and 1)",
     "", ""); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceModel::compute_temp_pa(IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = result.set_name("temp_pa"); CHKERRQ(ierr);
  ierr = result.set_attrs(
     "diagnostic",
     "pressure-adjusted ice temperature (degrees above pressure-melting point)",
     "deg_C", ""); CHKERRQ(ierr);

  PetscScalar *Tpaij, *Enthij; // columns of these values
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&Tpaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
	const PetscScalar p = EC->getPressureFromDepth(depth);
        ierr = EC->getPATemp(Enthij[k],p,Tpaij[k]);
          CHKERRQ(ierr);
	  if (config.get_flag("do_cold_ice_methods")) {
	    if ( EC->isTemperate(Enthij[k],p) && (vH(i,j) > 0)) {
	      Tpaij[k] = config.get("water_melting_temperature");
	    }
	  }
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // make deg C:
  ierr = result.shift(-config.get("water_melting_temperature")); CHKERRQ(ierr);

  // communication not done; we allow global IceModelVec3s as result
  return 0;
}

//! Computes the subglacial (basal) water pressure
/*!
  \f[p_w = \alpha\, \frac{w}{w_{\text{max}}}\, \rho\, g\, H,\f]
  where 

  - \f$\alpha\f$ is the till pore water fraction (till_pw_fraction),
  - \f$w\f$ is the effective thickness of subglacial melt water (bwat)
  - \f$w_{\text{max}}\f$ is the maximum allowed value for \f$w\f$ (hmelt_max),
  - \f$\rho\f$ is the ice density (ice_density)
  - \f$H\f$ is the ice thickness (thk)

Result is set to invalid (_FillValue) where the ice is floating, there being
no meaning to the above calculation.
 */
PetscErrorCode IceModel::compute_bwp(IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscScalar
    alpha     = config.get("till_pw_fraction"),
    wmax      = config.get("hmelt_max"),
    fillval   = -0.01;

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHmelt.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0.0) {
        result(i,j) = getBasalWaterPressure(
                        vH(i,j), vHmelt(i,j), vbmr(i,j), alpha, wmax);
      } else { // put negative value below valid range
        result(i,j) = fillval;
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = vMask.fill_where_floating(result, fillval); CHKERRQ(ierr);

  ierr = result.set_name("bwp"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "subglacial (pore) water pressure",
			  "Pa", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min",0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue",fillval); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes the rate of change of ice surface elevation as a sum of the
//! bedrock uplift rate and the thickness rate of change.
PetscErrorCode IceModel::compute_dhdt(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.copy_from(vdHdt); CHKERRQ(ierr); // result = dHdt
  ierr = result.mask_by(vH,0.0); CHKERRQ(ierr);	// set _FillValue areas to 0.0
  ierr = result.add(1.0, vuplift); CHKERRQ(ierr); // result += dbdt = uplift

  ierr = result.set_name("dhdt"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "rate of change of surface elevation",
			  "m s-1", ""); CHKERRQ(ierr);
  ierr = result.set_glaciological_units("m year-1");
  result.write_in_glaciological_units = true;

  return 0;
}

//! Computes the thickness of the basal layer of the temperate ice.
/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
PetscErrorCode IceModel::compute_tempicethk_basal(IceModelVec2S& result) {
  PetscErrorCode ierr;
  PetscScalar *Enth, fill_value = -0.01;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (vH(i,j) < 0.1)
        continue;

      ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
      PetscReal pressure;
      PetscInt ks = grid.kBelowHeight(vH(i,j)),
        k = 0;

      while (k <= ks) {
        pressure = EC->getPressureFromDepth(vH(i,j) - grid.zlevels[k]);

        if (EC->isTemperate(Enth[k],pressure))
          k++;
        else
          break;
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        result(i,j) = 0;
        continue;
      }
      
      // the whole column is temperate (except, possibly, some ice between
      // zlevels[ks] and the total thickness; we ignore it)
      if (k == ks + 1) {
        result(i,j) = grid.zlevels[ks];
        continue;
      }

      PetscReal 
        pressure_0 = EC->getPressureFromDepth(vH(i,j) - grid.zlevels[k-1]),
        dz = grid.zlevels[k] - grid.zlevels[k-1],
        slope1 = (Enth[k] - Enth[k-1]) / dz,
        slope2 = (EC->getEnthalpyCTS(pressure) - EC->getEnthalpyCTS(pressure_0)) / dz;
      
      if (slope1 != slope2) {
        result(i,j) = grid.zlevels[k-1] +
          (EC->getEnthalpyCTS(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        result(i,j) = PetscMax(result(i,j), grid.zlevels[k-1]);
        result(i,j) = PetscMin(result(i,j), grid.zlevels[k]);
      } else {
        SETERRQ4(1, "This should never happen: (i=%d, j=%d, k=%d, ks=%d)\n",
                    i, j, k, ks);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempicethk_basal"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "thickness of the basal layer of temperate ice",
			  "m", ""); CHKERRQ(ierr);

  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes thickness of the temperate ice layer (if there is any).
PetscErrorCode IceModel::compute_tempicethk(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar *Enth;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0.) {
	ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
	PetscScalar tithk = 0.;
	const PetscInt ks = grid.kBelowHeight(vH(i,j));
        
	for (PetscInt k=0; k<ks; ++k) {
          PetscReal pressure = EC->getPressureFromDepth(vH(i,j) - grid.zlevels[k]);

	  if (EC->isTemperate(Enth[k], pressure)) {
	    tithk += grid.zlevels[k+1] - grid.zlevels[k];
	  }
	}

	if (EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j) - grid.zlevels[ks]))) {
	  tithk += vH(i,j) - grid.zlevels[ks];
	}

	result(i,j) = tithk;
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempicethk"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "temperate ice thickness (total column content)",
			  "m", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice enthalpy at the base of ice.
PetscErrorCode IceModel::compute_enthalpybase(IceModelVec2S &result) {
  PetscErrorCode ierr;

  // put basal ice temperature in vWork2d[0]
  ierr = Enth3.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("enthalpybase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice enthalpy at the base of ice",
			  "J kg-1", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice temperature at the base of ice.
PetscErrorCode IceModel::compute_tempbase(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = compute_temp(vWork3d); CHKERRQ(ierr);

  // put basal ice temperature in vWork2d[0]
  ierr = vWork3d.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("tempbase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice temperature at the base of ice",
			  "Kelvin", ""); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes pressure-adjusted ice temperature at the base of ice.
PetscErrorCode IceModel::compute_temppabase(IceModelVec3 &hasPATemp,
                                            IceModelVec2S &result) {
  PetscErrorCode ierr;

  // put basal pressure-adjusted ice temperature in 2d result
  ierr = hasPATemp.getHorSlice(result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result.set_name("temppabase"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", 
                          "pressure-adjusted ice temperature at the base of ice",
			  "degrees Celsius", ""); CHKERRQ(ierr);

  PetscScalar fill_value = 0.01;
  ierr = result.mask_by(vH, fill_value); CHKERRQ(ierr);
  ierr = result.set_attr("valid_max", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);

  return 0;
}

//! Computes ice temperature at the 1 m below the surface.
PetscErrorCode IceModel::compute_tempsurf(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01;

  // compute levels corresponding to 1 m below the ice surface:

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      result(i,j) = PetscMax(vH(i,j) - 1.0, 0.0);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = compute_temp(vWork3d); CHKERRQ(ierr);

  ierr = vWork3d.getSurfaceValues(result, result); CHKERRQ(ierr);  // z=H-1 slice

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) <= 1.0)
	result(i,j) = fill_value;
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice temperature at 1m below the ice surface",
			  "Kelvin", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);
  return 0;
}

//! Computes ice enthalpy at the 1 m below the surface.
PetscErrorCode IceModel::compute_enthalpysurf(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscScalar fill_value = -0.01;

  // compute levels corresponding to 1 m below the ice surface:

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      result(i,j) = PetscMax(vH(i,j) - 1.0, 0.0);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = Enth3.getSurfaceValues(result, result); CHKERRQ(ierr);  // z=0 slice

  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) <= 1.0)
	result(i,j) = fill_value;
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("tempsurf"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic", "ice enthalpy at 1m below the ice surface",
			  "J kg-1", ""); CHKERRQ(ierr);
  ierr = result.set_attr("_FillValue", fill_value); CHKERRQ(ierr);
  return 0;
}

//! \brief Computes a diagnostic quantity given by \c name and returns a
//! pointer to a pre-allocated work vector containing it.
/*! 
For 2D quantities, result will point to vWork2d[0].  For 3D -- to vWork3d.
Depending on the quantity requested, vWork2d[1] might get used as temporary
storage.
 */
PetscErrorCode IceModel::compute_by_name(string name, IceModelVec* &result) {
  PetscErrorCode ierr;

  result = NULL;		// if clauses can override this

  if (name == "bwp") {
    ierr = compute_bwp(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "cts") {
    ierr = compute_cts(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "dhdt") {
    ierr = compute_dhdt(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "enthalpybase") {
    ierr = compute_enthalpybase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "enthalpysurf") {
    ierr = compute_enthalpysurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "liqfrac") {
    ierr = compute_liqfrac(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "temp") {
    ierr = compute_temp(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "temp_pa") {
    ierr = compute_temp_pa(vWork3d); CHKERRQ(ierr);
    result = &vWork3d;
    return 0;
  }

  if (name == "tempsurf") {
    ierr = compute_tempsurf(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "tempbase") {
    ierr = compute_tempbase(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "temppabase") {
    ierr = compute_temp_pa(vWork3d); CHKERRQ(ierr);
    ierr = compute_temppabase(vWork3d,vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "tempicethk") {
    ierr = compute_tempicethk(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "tempicethk_basal") {
    ierr = compute_tempicethk_basal(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "proc_ice_area") {
    ierr = compute_proc_ice_area(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  if (name == "bueler_brown_f") {
    ierr = compute_bueler_brown_f(vWork2d[0]); CHKERRQ(ierr);
    result = &vWork2d[0];
    return 0;
  }

  return 0;
}


//! Computes the ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0)
        volume += vH(i,j) * cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Computes the temperate ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume_temperate(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  
  PetscScalar *Enth;  // do NOT delete this pointer: space returned by
                      //   getInternalColumn() is allocated already
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
        ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        for (PetscInt k=0; k<ks; ++k) {	  
          if (EC->isTemperate(Enth[k],EC->getPressureFromDepth(vH(i,j)))) {
            volume += (grid.zlevels[k+1] - grid.zlevels[k]) * cell_area(i,j);
          }
        }
        if (EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j)))) {
          volume += (vH(i,j) - grid.zlevels[ks]) * cell_area(i,j);
        }
      }
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes the cold ice volume, in m^3.
PetscErrorCode IceModel::compute_ice_volume_cold(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     volume=0.0;
  
  PetscScalar *Enth;  // do NOT delete this pointer: space returned by
                      //   getInternalColumn() is allocated already
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
        ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        for (PetscInt k=0; k<ks; ++k) {	  
          if (!EC->isTemperate(Enth[k],EC->getPressureFromDepth(vH(i,j)))) {
            volume += (grid.zlevels[k+1] - grid.zlevels[k]) * cell_area(i,j);
          }
        }
        if (!EC->isTemperate(Enth[ks],EC->getPressureFromDepth(vH(i,j)))) {
          volume += (vH(i,j) - grid.zlevels[ks]) * cell_area(i,j);
        }
      }
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&volume, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0)
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes temperate ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_temperate(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && (EC->isTemperate(Enthbase[i][j],EC->getPressureFromDepth(vH(i,j)))) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes cold ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_cold(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && (!EC->isTemperate(Enthbase[i][j],EC->getPressureFromDepth(vH(i,j)))) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes grounded ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_grounded(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && vMask.is_grounded(i,j) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Computes floating ice area, in m^2.
PetscErrorCode IceModel::compute_ice_area_floating(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar     area=0.0;
  
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ( (vH(i,j) > 0) && vMask.is_floating(i,j) )
        area += cell_area(i,j);
    }
  }  
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&area, &result, grid.com); CHKERRQ(ierr);
  return 0;
}


//! Computes the total ice enthalpy in J.
/*!
Units of the specific enthalpy field \f$E=\f$(IceModelVec3::Enth3) are J kg-1.  We integrate
\f$E(t,x,y,z)\f$ over the entire ice fluid region \f$\Omega(t)\f$, multiplying
by the density to get units of energy:
   \f[ E_{\text{total}}(t) = \int_{\Omega(t)} E(t,x,y,z) \rho_i \,dx\,dy\,dz. \f]
 */
PetscErrorCode IceModel::compute_ice_enthalpy(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar enthalpysum = 0.0;

  PetscScalar *Enth;  // do NOT delete this pointer: space returned by
                      //   getInternalColumn() is allocated already
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
        ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
        for (PetscInt k=0; k<ks; ++k) {
          enthalpysum += Enth[k] * (grid.zlevels[k+1] - grid.zlevels[k]);
        }
        enthalpysum += Enth[ks] * (vH(i,j) - grid.zlevels[ks]);
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);  
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  
  enthalpysum *= config.get("ice_density") * (grid.dx * grid.dy);
  
  ierr = PetscGlobalSum(&enthalpysum, &result, grid.com); CHKERRQ(ierr);
  return 0;
}

//! Compute a scalar diagnostic quantity by name.
PetscErrorCode IceModel::compute_by_name(string name, PetscScalar &result) {
  PetscErrorCode ierr, errcode = 1;

  if (name == "ivol") {
    errcode = 0;
    ierr = compute_ice_volume(result); CHKERRQ(ierr);
  }

  if (name == "ivoltemp") {
    errcode = 0;
    ierr = compute_ice_volume_temperate(result); CHKERRQ(ierr);
  }

  if (name == "ivoltempf") {
    errcode = 0;
    PetscScalar ivol;
    ierr = compute_ice_volume(ivol); CHKERRQ(ierr);
    ierr = compute_ice_volume_temperate(result); CHKERRQ(ierr);
    result /= ivol;
  }

  if (name == "ivolcold") {
    errcode = 0;
    ierr = compute_ice_volume_cold(result); CHKERRQ(ierr);
  }

  if (name == "ivolcoldf") {
    errcode = 0;
    PetscScalar ivol;
    ierr = compute_ice_volume(ivol); CHKERRQ(ierr);
    ierr = compute_ice_volume_cold(result); CHKERRQ(ierr);
    result /= ivol;
  }

  if (name == "imass") {
    errcode = 0;
    PetscScalar ice_density = config.get("ice_density");
    ierr = compute_ice_volume(result); CHKERRQ(ierr);
    result *= ice_density;
  }

  if (name == "iarea") {
    errcode = 0;
    ierr = compute_ice_area(result); CHKERRQ(ierr);
  }

  if (name == "iareatemp") {
    errcode = 0;
    ierr = compute_ice_area_temperate(result); CHKERRQ(ierr);
  }

  if (name == "iareatempf") {
    errcode = 0;
    PetscScalar iarea;
    ierr = compute_ice_area(iarea); CHKERRQ(ierr);
    ierr = compute_ice_area_temperate(result); CHKERRQ(ierr);
    result /= iarea;
  }

  if (name == "iareacold") {
    errcode = 0;
    ierr = compute_ice_area_cold(result); CHKERRQ(ierr);
  }

  if (name == "iareacoldf") {
    errcode = 0;
    PetscScalar iarea;
    ierr = compute_ice_area(iarea); CHKERRQ(ierr);
    ierr = compute_ice_area_cold(result); CHKERRQ(ierr);
    result /= iarea;
  }

  if (name == "iareag") {
    errcode = 0;
    ierr = compute_ice_area_grounded(result); CHKERRQ(ierr);
  }

  if (name == "iareaf") {
    errcode = 0;
    ierr = compute_ice_area_floating(result); CHKERRQ(ierr);
  }

  if (name == "dt") {
    errcode = 0;
    result = dt;
  }

  if (name == "divoldt") {
    errcode = 0;
    result = dvoldt;
  }

  if (name == "dimassdt") {
    errcode = 0;
    PetscScalar ice_density = config.get("ice_density");
    result = dvoldt * ice_density;
  }

  if (name == "ienthalpy") {
    errcode = 0;
    ierr = compute_ice_enthalpy(result); CHKERRQ(ierr);
  }

  if (name == "basal_ice_flux") {
    errcode = 0;
    result = total_basal_ice_flux;
  }

  if (name == "surface_ice_flux") {
    errcode = 0;
    result = total_surface_ice_flux;
  }

  if (name == "sub_shelf_ice_flux") {
    errcode = 0;
    result = total_sub_shelf_ice_flux;
  }

  if (name == "nonneg_rule_flux") {
    errcode = 0;
    result = nonneg_rule_flux;
  }

  if (name == "ocean_kill_flux") {
    errcode = 0;
    result = ocean_kill_flux;
  }

  if (name == "float_kill_flux") {
    errcode = 0;
    result = float_kill_flux;
  }

  return errcode;
}


/*! Computes total ice fluxes in kg s-1 at 3 interfaces:

  \li the ice-atmosphere interface: gets surface mass balance rate from
      PISMSurfaceModel *surface,
  \li the ice-ocean interface at the bottom of ice shelves: gets ocean-imposed
      basal melt rate from PISMOceanModel *ocean, and
  \li the ice-bedrock interface: gets basal melt rate from IceModelVec2S vbmr.

A unit-conversion occurs for all three quantities, from ice-equivalent m s-1
to kg s-1.  The sign convention about these fluxes is that positve flux means
ice is being \e added to the ice fluid volume at that interface.

These quantities should be understood as <i>instantaneous at the beginning of
the time-step.</i>  Multiplying by dt will \b not necessarily give the
corresponding change from the beginning to the end of the time-step.

FIXME:  The calving rate can be computed by post-processing:
divoldt = surface_ice_flux * iarea + basal_ice_flux * iareag + sub_shelf_ice_flux * iareaf + calving_flux_vol_rate
 */
PetscErrorCode IceModel::ice_mass_bookkeeping() {
  PetscErrorCode ierr;

  bool include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity");

  // note acab and shelfbmassflux are IceModelVec2S owned by IceModel
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_mass_flux(grid.year, dt / secpera, acab);
    CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: surface == PETSC_NULL"); }

  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(grid.year, dt / secpera, shelfbmassflux);
    CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: ocean == PETSC_NULL"); }

  PetscScalar my_total_surface_ice_flux = 0.0, my_total_basal_ice_flux = 0.0,
    my_total_sub_shelf_ice_flux = 0.0;

  ierr = acab.begin_access(); CHKERRQ(ierr);
  ierr = cell_area.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // ignore ice-free cells:
      if (vH(i,j) <= 0.0)
        continue;

      my_total_surface_ice_flux += acab(i,j) * cell_area(i,j); // note the "+="!

      if ((vMask.value(i,j) == MASK_FLOATING) && include_bmr_in_continuity) {
        // note: we are deliberately *not* including fluxes in
        //   MASK_ICE_FREE_OCEAN and MASK_OCEAN_AT_TIME_0 areas
        my_total_sub_shelf_ice_flux -= shelfbmassflux(i,j) * cell_area(i,j); // note the "-="!
      }

      if (vMask.is_grounded(i,j) && include_bmr_in_continuity) {
        my_total_basal_ice_flux -= vbmr(i,j) * cell_area(i,j); // note the "-="!
      }
    }	// j
  } // i

  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = cell_area.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  PetscScalar ice_density = config.get("ice_density");
  my_total_surface_ice_flux     *= ice_density;
  my_total_sub_shelf_ice_flux   *= ice_density;
  my_total_basal_ice_flux       *= ice_density;

  ierr = PetscGlobalSum(&my_total_surface_ice_flux,   &total_surface_ice_flux,   grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&my_total_sub_shelf_ice_flux, &total_sub_shelf_ice_flux, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&my_total_basal_ice_flux,     &total_basal_ice_flux,     grid.com); CHKERRQ(ierr);

  return 0;
}

