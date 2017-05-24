#ifndef MJOLNIR_LOCAL_INTEARACTION_BASE
#define MJOLNIR_LOCAL_INTEARACTION_BASE
#include "System.hpp"
#include <array>

namespace mjolnir
{

template<typename traitsT>
class LocalInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename system_type::particle_type   particle_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void
    calc_force(system_type&) const = 0;

    virtual real_type
    calc_energy(const system_type&) const = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
