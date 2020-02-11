#ifndef __LandEvo_hh
#define __LandEvo_hh


#include "BedDef.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace bed {

class LandEvo : public BedDef {
public:
  LandEvo(IceGrid::ConstPtr g);
  ~LandEvo();

  void init(const InputOptions &opts);
  void update(
    const IceModelVec3 &u3,
    const IceModelVec3 &v3,
    const IceModelVec2CellType &mask,
    double dt);

protected:
  IceModelVec2S m_topg;
  IceModelVec2S m_prescribed_uplift;
  bool do_erosion, do_prescribed_uplift;

  void init_impl(const InputOptions &opts);
  void update_impl(const IceModelVec2S &ice_thickness,
                   const IceModelVec2S &sea_level_elevation,
                   double t, double dt);
  void update_erosion(
    const IceModelVec3 &u3,
    const IceModelVec3 &v3,
    const IceModelVec2CellType &mask,
    double dt);
  void update_lem(
    const IceModelVec3 &u3,
    const IceModelVec3 &v3,
    const IceModelVec2CellType &mask,
    double dt);
  void update_prescribed_uplift(double dt);
  void update_with_thickness_impl(
    const IceModelVec2S &ice_thickness,
    double my_t, double my_dt);
};

}
}

#endif  //__BEDLEM_HH