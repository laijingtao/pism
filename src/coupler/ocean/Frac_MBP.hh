/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef _FRAC_MBP_H_
#define _FRAC_MBP_H_

#include "pism/coupler/util/PScalarForcing.hh"
#include "pism/coupler/OceanModel.hh"

namespace pism {
namespace ocean {

/**
 * Scalar melange back-pressure fraction forcing.
 * 
 */
class Frac_MBP : public PScalarForcing<OceanModel,OceanModel>
{
public:
  Frac_MBP(IceGrid::ConstPtr g, OceanModel* in);
  virtual ~Frac_MBP();

protected:
  virtual void init_impl();

  virtual void update_impl(double t, double dt);

  const IceModelVec2S& melange_back_pressure_fraction_impl() const;
private:
  typedef PScalarForcing<OceanModel,OceanModel> super;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _FRAC_MBP_H_ */
