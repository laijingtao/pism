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
  void update_lem(
    const IceModelVec2S &ice_thickness,
    const IceModelVec3 &u3,
    const IceModelVec3 &v3,
    const IceModelVec2CellType &mask,
    double dt);

protected:
  bool do_erosion = false, do_prescribed_uplift = false, do_fluvial = false;
  IceModelVec2S m_prescribed_uplift;

  void init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                 const IceModelVec2S &sea_level_elevation);
  void update_impl(const IceModelVec2S &ice_thickness,
                   const IceModelVec2S &sea_level_elevation,
                   double t, double dt);
  void update_erosion(
    const IceModelVec2S &ice_thickness,
    const IceModelVec3 &u3,
    const IceModelVec3 &v3,
    const IceModelVec2CellType &mask,
    double dt);
  void update_prescribed_uplift(double dt);
  void compute_erosion_threshold(const std::string &stabilizing_method,
                                 const IceModelVec2S &bed_elevation,
                                 const IceModelVec2S &ice_thickness,
                                 IceModelVec2S &results);
  void update_fluvial_erosion(
    const IceModelVec2S &ice_thickness,
    const IceModelVec2CellType &mask,
    double dt);
  void fluvial_route_flow(
    const IceModelVec2S &bed_elevation,
    const IceModelVec2S &ice_thickness,
    std::vector<std::vector<double>> &drainage_area,
    std::vector<std::vector<double>> &steepest_slope);
};

struct node_coord {
  int i, j;
};


void add_to_stack(int i, int j,
    const std::vector<std::vector<node_coord>> &receiver,
    std::vector<node_coord> &ordered_nodes,
    std::vector<std::vector<bool>> &in_list);

}
}

#endif  //__BEDLEM_HH