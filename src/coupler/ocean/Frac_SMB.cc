// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <gsl/gsl_math.h>

#include "Frac_SMB.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace ocean {

Frac_SMB::Frac_SMB(IceGrid::ConstPtr g, OceanModel* in)
  : PScalarForcing<OceanModel,OceanModel>(g, in) {

  m_option_prefix = "-ocean_frac_mass_flux";
  m_offset_name   = "frac_mass_flux";

  m_offset.reset(new Timeseries(*m_grid, m_offset_name, m_config->get_string("time.dimension_name")));

  m_offset->variable().set_string("units", "1");
  m_offset->dimension().set_string("units", m_grid->ctx()->time()->units_string());
  m_offset->variable().set_string("long_name", "ice-shelf-base mass flux factor");

  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);
}

Frac_SMB::~Frac_SMB() {
  // empty
}

void Frac_SMB::init_impl() {

  m_input_model->init();

  m_log->message(2,
             "* Initializing ice shelf base mass flux forcing using scalar offsets...\n");

  init_internal();
}

void Frac_SMB::update_impl(double t, double dt) {
  super::update_impl(t, dt);

  scale_data(*m_shelf_base_mass_flux);
}

const IceModelVec2S& Frac_SMB::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}

} // end of namespace ocean
} // end of namespace pism
