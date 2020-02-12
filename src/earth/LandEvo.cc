// Copyright (C) 2010, 2011, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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

#include "BedDef.hh"
#include "LandEvo.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Time.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Vars.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace bed {

LandEvo::LandEvo(IceGrid::ConstPtr g)
  : BedDef(g) {

  m_prescribed_uplift.create(m_grid, "prescribed_uplift", WITHOUT_GHOSTS);
  m_prescribed_uplift.set_attrs("model_state", "prescribed bedrock uplift rate",
                           "m year-1", "prescribed_uplift_rate");

  do_erosion = m_config->get_boolean("bed_deformation.erosion.enabled");
  do_prescribed_uplift = m_config->get_boolean("bed_deformation.prescribed_uplift.enabled");
}

LandEvo::~LandEvo() {
  // empty
}

void LandEvo::init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                        const IceModelVec2S &sea_level_elevation) {

  // ice_thickness and sea_level_elevation are not used (for now)
  (void) ice_thickness;
  (void) sea_level_elevation;

  land_evo_model_enabled = true;

  m_log->message(2,
             "* Initializing the landscape evolution model (bed deformation model)...\n");

  if (do_erosion) {
    m_log->message(2, "   - glacial erosion is enabled \n");
  }

  if (do_prescribed_uplift) {
    m_log->message(2, "   - prescribed uplift is enabled \n");
  }

  BedDef::init_impl(opts, ice_thickness, sea_level_elevation);

  m_uplift.set(0.0);

  if (do_prescribed_uplift) {
    std::string prescribed_uplift_file = m_config->get_string("bed_deformation.prescribed_uplift_file");
    if (prescribed_uplift_file.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "missing prescribed_uplift_file");
    }
    m_log->message(2,
                   "    reading bed uplift from %s ... \n",
                   prescribed_uplift_file.c_str());
    m_prescribed_uplift.regrid(prescribed_uplift_file, CRITICAL);
  }
}

void LandEvo::update_impl(const IceModelVec2S &ice_thickness,
                       const IceModelVec2S &sea_level_elevation,
                       double t, double dt) {

  // have to have this func doing nothing because it is a pure virtual func in BedDef.hh
  (void) ice_thickness;
  (void) sea_level_elevation;
  (void) t;
  (void) dt;
}

void LandEvo::update_lem(
    const IceModelVec3 &u3,
    const IceModelVec3 &v3,
    const IceModelVec2CellType &mask,
    double dt) {

  if (do_erosion) {
    this->update_erosion(u3, v3, mask, dt);
  }

  if (do_prescribed_uplift) {
    this->update_prescribed_uplift(dt);
  }
}

void LandEvo::update_erosion(const IceModelVec3 &u3,
                            const IceModelVec3 &v3,
                            const IceModelVec2CellType &mask,
                            double dt) {
  IceModelVec2S tmp_u, tmp_v;
  tmp_u.create(m_grid, "tmp_u", WITHOUT_GHOSTS);
  tmp_v.create(m_grid, "tmp_v", WITHOUT_GHOSTS);

  u3.getHorSlice(tmp_u, 0.0);
  v3.getHorSlice(tmp_v, 0.0);

  IceModelVec2S::Ptr sliding_mag_ptr(new IceModelVec2S(m_grid, "sliding_mag", WITHOUT_GHOSTS));
  sliding_mag_ptr->set_to_magnitude(tmp_u, tmp_v);
  IceModelVec2S &sliding_mag = *sliding_mag_ptr;

  double k_g = m_config->get_double("bed_deformation.erosion.k_g");
  double l = m_config->get_double("bed_deformation.erosion.exponent");

  IceModelVec::AccessList list{&m_topg, &sliding_mag, &mask};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if (mask.grounded_ice(i, j)) {
        if (l == 1.0) {
          m_topg(i, j) = m_topg(i, j) - dt / 31557600.0 * (k_g * sliding_mag(i, j) * 31557600.0);
        } else {
          m_topg(i, j) = m_topg(i, j) - dt / 31557600.0 * (k_g * pow(sliding_mag(i, j) * 31557600.0, l));
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_topg.inc_state_counter();
}

void LandEvo::update_prescribed_uplift(double dt) {
  IceModelVec::AccessList list{&m_topg, &m_prescribed_uplift};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      m_topg(i, j) = m_topg(i, j) + m_prescribed_uplift(i, j) * dt / 31557600.0;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_topg.inc_state_counter();
}


} // end of namespace bed
} // end of namespace pism
