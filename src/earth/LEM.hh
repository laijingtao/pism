#ifndef __BEDLEM_HH
#define __BEDLEM_HH


#include "BedDef.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace bed {

class LEMInputs {
public:
    LEMInputs();
    LEMInputs(const IceModelVec3 *u3,
              const IceModelVec3 *v3,
              const IceModelVec2CellType *mask,
              double dt);

    const IceModelVec3 *u3, *v3;
    const IceModelVec2CellType *mask;
    double dt;
};

class LandscapeEvolution : public Component {
public:
    LandscapeEvolution(IceGrid::ConstPtr g);
    virtual ~LandscapeEvolution();

    void update_lem(LEMInputs inputs);

protected:
    IceModelVec2S m_fixed_uplift;

    void init_impl(const InputOptions &opts);
    void update_erosion(
        const IceModelVec3 &u3,
        const IceModelVec3 &v3,
        const IceModelVec2CellType &mask,
        double dt);
    void update_fixed_uplift(double dt);
    void update_with_thickness_impl(
        const IceModelVec2S &ice_thickness,
        double my_t, double my_dt);
};

}
}

#endif  //__BEDLEM_HH