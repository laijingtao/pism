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
#include "LEM.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Time.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Vars.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace bed {

LandscapeEvolution::LandscapeEvolution(IceGrid::ConstPtr g)
  : BedDef(g) {
  m_fixed_uplift.create(m_grid, "fixed_uplift", WITHOUT_GHOSTS);
  m_fixed_uplift.set_attrs("model_state", "fixed bedrock uplift rate",
                           "m year-1", "fixed_uplift_rate");
}

void LandscapeEvolution::init_impl(const InputOptions &opts) {

  m_log->message(2,
             "* Initializing the landscape evolution model...\n");

  BedDef::init_impl(opts);

  bool do_fixed_uplift = m_config->get_boolean("bed_deformation.fixed_uplift.enabled");
  if (do_fixed_uplift) {
    std::string fixed_uplift_file = m_config->get_string("bed_deformation.fixed_uplift_file");
    if (fixed_uplift_file.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "missing fixed_uplift_file");
    }
    m_log->message(2,
                   "    reading bed uplift from %s ... \n",
                   fixed_uplift_file.c_str());
    m_fixed_uplift.regrid(fixed_uplift_file, CRITICAL);
  }

}

void LandscapeEvolution::update_lem(LEMInputs inputs) {
  bool do_erosion = m_config->get_boolean("bed_deformation.erosion.enabled");
  if (do_erosion) {
    this->update_erosion(*(inputs.u3), *(inputs.v3), *(inputs.mask), inputs.dt);
  }

  bool do_fixed_uplift = m_config->get_boolean("bed_deformation.fixed_uplift.enabled");
  if (do_fixed_uplift) {
    this->update_fixed_uplift(inputs.dt)
  }
}

void LandscapeEvolution::update_erosion(const IceModelVec3 &u3,
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
  IceModelVec2S sliding_mag = *sliding_mag_ptr;

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

void LandscapeEvolution::update_fixed_uplift(double dt) {
  IceModelVec::AccessList list{&m_topg, &m_fixed_uplift};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      m_topg(i, j) = m_topg(i, j) + m_fixed_uplift(i, j) * dt / 31557600.0;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_topg.inc_state_counter();
}

void LandscapeEvolution::update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                  double my_t, double my_dt) {
  (void) ice_thickness;
  (void) t;
  (void) dt;
  // This model does not update isostatic change.
}

} // end of namespace bed
} // end of namespace pism
