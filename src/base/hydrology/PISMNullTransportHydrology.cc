// Copyright (C) 2012-2013 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "PISMHydrology.hh"
#include "hydrology_diagnostics.hh"


PetscErrorCode PISMNullTransportHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the null-transport (till only) subglacial hydrology model ...\n"); CHKERRQ(ierr);
  ierr = PISMHydrology::init(vars); CHKERRQ(ierr);
  return 0;
}


//! Set the transportable subglacial water thickness to zero; there is no tranport.
PetscErrorCode PISMNullTransportHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}


//! Returns the (trivial) overburden pressure as the pressure of the non-existent transportable water, because this is the least harmful output if this is misused.
PetscErrorCode PISMNullTransportHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr = overburden_pressure(result); CHKERRQ(ierr);
  return 0;
}


//! Update the till water thickness by simply integrating the melt input.  There is no tranportable water thickness variable and no interaction with it.
/*!
Does an explicit (Euler) step of the integration
  \f[ \frac{\partial W_{til}}{\partial t} = \min \left\{\frac{m}{\rho_w} - C\f]
where \f$C=\f$`hydrology_tillwat_decay_rate_null`.  Enforces bounds
\f$0 \le W_{til} \le W_{til}^{max}\f$ where the upper bound is
`hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.
 */
PetscErrorCode PISMNullTransportHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;

  ierr = get_input_rate(icet,icedt,total_input); CHKERRQ(ierr);

  const PetscReal Wtilmax  = config.get("hydrology_tillwat_max"),
                  C        = config.get("hydrology_tillwat_decay_rate_null");
  if (Wtilmax < 0.0) {
    PetscPrintf(grid.com,
       "PISMNullTransportHydrology ERROR: hydrology_tillwat_max is negative\n"
       "            this is not allowed ... ENDING ... \n\n");
    PISMEnd();
  }

  ierr = Wtil.begin_access(); CHKERRQ(ierr);
  ierr = total_input.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Wtil(i,j) += icedt * (total_input(i,j) - C);
      Wtil(i,j) = PetscMin(PetscMax(0.0, Wtil(i,j)), Wtilmax);
    }
  }
  ierr = Wtil.end_access(); CHKERRQ(ierr);
  ierr = total_input.end_access(); CHKERRQ(ierr);
  return 0;
}

