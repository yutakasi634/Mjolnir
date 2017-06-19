#ifndef MJOLNIR_GLOBAL_BURIEDNESS_INTEARACTION
#define MJOLNIR_GLOBAL_BURIEDNESS_INTEARACTION
#include "GlobalInteractionBase.hpp"
#include <queue>

namespace mjolnir
{

/* This Interaction class is especially for Hydrophobic Interaction developed *
 * by Tanaka et al. (2015), which has a similar functional form to that was   *
 * used in Fujitsuka et al. (2004).                                           *
 * It requires Potential class to have `rmin(i, j)`, `rmax(i, j)`m,           *
 * `buriedness(i, j, dist_ij)` and its derivative, `dburiedness(i, j, dij)`   *
 * in addition to normal `calc_force` and `calc_energy`.                      *
 * So this is an ad hoc Interaction class that is possibly changed later.     */
template<typename traitsT, typename potentialT, typename partitionT>
class BuriednessInteraction final : public GlobalInteractionBase<traitsT>
{
  public:
    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef partitionT partition_type;
    typedef GlobalInteractionBase<traitsT> base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::particle_type   particle_type;
    typedef typename base_type::boundary_type   boundary_type;

  public:
    BuriednessInteraction()  = default;
    ~BuriednessInteraction() = default;

    BuriednessInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    void initialize(const system_type& sys, const real_type dt) override
    {
        tmpforce_.reserve(20);
        this->partition_.initialize(sys);
        this->partition_.update(sys, dt);
    }

    void      calc_force (system_type&)             override;
    real_type calc_energy(const system_type&) const override;

  private:

    potential_type potential_;
    partition_type partition_;
    std::vector<std::pair<std::size_t, coordinate_type>> tmpforce_;
};

template<typename traitsT, typename potT, typename spaceT>
void BuriednessInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys)
{
    potential_.update(sys);
    partition_.update(sys);

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        real_type rho(0.), duhp(0.);
        tmpforce_.emplace_back(i, coordinate_type(0., 0., 0.));

        for(auto j : this->partition_.partners(i))
        {
            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - sys[i].position);
            const real_type l = length(rij);

            rho += potential_.buriedness(i, j, l);
            const coordinate_type f = potential_.dburiedness(i, j, l) * rij;

            tmpforce_.front() += f;
            tmpforce_.emplace_back(j, -1 * f);
        }

        const real_type coef = potential_.derivative(i, rho);
        for(const auto& tf : tmpforce_)
        {
            sys[tf.first].force += coef * tf.second;
        }
        tmpforce_.clear();
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename BuriednessInteraction<traitsT, potT, spaceT>::real_type
BuriednessInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        real_type rho = 0.;
        for(auto j : this->partition_.partners(i))
        {
            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - sys[i].position);
            rho += potential_.buriedness(i, j, length(rij));
        }
        e += potential_.potential(i, rho);
    }
    return e;
}


} // mjolnir
#endif//MJOLNIR_GLOBAL_BURIEDNESS_INTEARACTION
