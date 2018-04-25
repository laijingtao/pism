// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Andy Aschwanden and Constantine Khroulev
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

#include "Avalanche.hh"

#include "pism/util/iceModelVec.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace surface {

///// Avalanche surface model.
Avalanche::Avalanche(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> in)
  : SurfaceModel(g) {

  m_mass_flux   = allocate_mass_flux(g);
}

void Avalanche::init_impl(const Geometry &geometry) {
  (void) geometry;

  bool limits_set = false;

  m_log->message(2,
                 "* Initializing the surface processes model Avalanche. Setting...\n");

}

MaxTimestep Avalanche::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("surface 'avalanche'");
}

void Avalanche::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  compute_mass_flux(geometry.ice_surface_elevation, *m_mass_flux);
}

const IceModelVec2S &Avalanche::mass_flux_impl() const {
  return *m_mass_flux;
}


void Avalanche::compute_mass_flux(const IceModelVec2S &surface, IceModelVec2S &result) const {
  double dabdz = -m_M_min/(m_z_ELA - m_z_M_min);
  double dacdz = m_M_max/(m_z_M_max - m_z_ELA);

  IceModelVec::AccessList list{&result, &surface};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double z = surface(i, j);

      if (z < m_z_M_min) {
        result(i, j) = m_M_limit_min;
      }
      else if ((z >= m_z_M_min) and (z < m_z_ELA)) {
        result(i, j) = dabdz * (z - m_z_ELA);
      }
      else if ((z >= m_z_ELA) and (z <= m_z_M_max)) {
        result(i, j) = dacdz * (z - m_z_ELA);
      }
      else if (z > m_z_M_max) {
        result(i, j) = m_M_limit_max;
      }
      else {
        throw RuntimeError(PISM_ERROR_LOCATION,
                           "Avalanche::compute_mass_flux: HOW DID I GET HERE?");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // convert from m second-1 ice equivalent to kg m-2 s-1:
  result.scale(m_config->get_double("constants.ice.density"));
}


} // end of namespace surface
} // end of namespace pism
