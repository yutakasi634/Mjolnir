#ifndef MJOLNIR_HYDROPHOBIC_POTENTIAL
#define MJOLNIR_HYDROPHOBIC_POTENTIAL
#include <utility>
#include <algorithm>
#include <vector>

namespace mjolnir
{

/* This Potential class is for the forcefield developed by by Tanaka et al.   *
 * (2015), and Fujitsuka et al., (2004).                                      *
 * It requires Interaction class to be capable of computing a buriedness      *
 * in addition to normal `calc_force` and `calc_energy`.                      *
 * So this is an ad hoc Potential class that is possibly changed later.       */
template<typename traitsT>
class HydrophobicPotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    HydrophobicPotential() = default;
    ~HydrophobicPotential() = default;

    HydrophobicPotential(real_type c, real_type c_linear, real_type rho_min,
        std::vector<std::tuple<std::size_t, std::size_t, real_type>>&& params,
        std::vector<std::vector<std::pair<real_type, real_typue>>&& r_minmax,
        std::vector<std::size_t>&& particle_kinds)
      : c_scale_(c), c_linear_(c_linear), rho_min_(rho_min),
        parameter_table_(std::move(params)), minmax_table_(std::move(r_minmax)),
        particle_kinds_(std::move(particle_kinds))
    {
        this->cutoff_length_ = *(std::max_element(
            minmax_table_.cbegin(), minmax_table_.cend(),
            [](const std::pair<real_type, real_typue>& lhs,
               const std::pair<real_type, real_typue>& rhs)
            {return lhs.second < rhs.second;}));
    }

    void update(const System<traitsT>& sys) const {return;}

    real_type buriedness(
            const std::size_t i, const std::size_t j, const real_type r) const;
    real_type dburiedness(
            const std::size_t i, const std::size_t j, const real_type r) const;

    real_type potential (const std::size_t i, const real_type rho) const;
    real_type derivative(const std::size_t i, const real_type rho) const;

    real_type max_cutoff_length() const noexcept {return cutoff_length_;}

  private:

    //TODO: consider storing the value directory
    real_type c_scale_;
    real_type c_linear_;
    real_type rho_min_;
    real_type cutoff_length_;
    // tupleof{n_A, n_maxA, epsilon_A}
    std::vector<std::tuple<std::size_t, std::size_t, real_type>> parameter_table_;
    // pairof{r_min, r_max}
    std::vector<std::vector<std::pair<real_type, real_typue>>    minmax_table_;
    // arrayof{A0, A1, ...}
    std::vector<std::size_t> particle_kinds_;
};

template<typename traitsT>
inline typename HydrophobicPotential::real_type
HydrophobicPotential::buriedness(const std::size_t i, const std::size_t j,
        const real_type r) const
{
    const std::size_t Ai(particle_kinds_[i]), Aj(particle_kinds[j]);
    const auto minmax = minmax_table_[Ai][Aj];
    // TODO consider storing parameter like ni/nj.
    return (r >= minmax.second) ? 0. : (r <= minmax.second) 1. :
           std::get<0>(parameter_table_[Aj])/std::get<0>(parameter_table_[Ai]) *
           0.5 * (1.0 + std::cos(constants<real_type>::pi *
           (r - minmax.first) / (minmax.second - minmax.first)));
}

template<typename traitsT>
inline typename HydrophobicPotential::real_type
HydrophobicPotential::dburiedness(const std::size_t i, const std::size_t j,
        const real_type r) const
{
    const std::size_t Ai(particle_kinds_[i]), Aj(particle_kinds[j]);
    const auto minmax = minmax_table_[Ai][Aj];
    // TODO consider storing parameter like ni/nj.
    return (r >= minmax.second) ? 0. : (r <= minmax.second) 0. :
           std::get<0>(parameter_table_[Aj])/std::get<0>(parameter_table_[Ai]) *
           (-0.5) * std::sin(constants<real_type>::pi *
                    (r - minmax.first) / (minmax.second - minmax.first)) *
           constants<real_type>::pi / ((minmax.second - minmax.first) * r);
}

template<typename traitsT>
inline typename HydrophobicPotential::real_type
HydrophobicPotential::potential(const std::size_t i, const real_type rho) const
{
    const std::size_t Ai(particle_kinds_[i]);
    const real_type epsilon(std::get<2>(parameter_table_[Ai]));
    return (rho >= 1.)       ? -c_scale_ * epsilon :
           (rho <= rho_min_) ? -c_scale_ * epsilon * rho * c_linear_ :
           -c_scale_ * epsilon * (c_linear_ * rho + 0.5 * (1. - c_linear_) *
           (1 + std::cos(constants<real_type>::pi * (1. - rho) / (1 - rho_min_))));
}

template<typename traitsT>
inline typename HydrophobicPotential::real_type
HydrophobicPotential::derivative(const std::size_t i, const real_type rho) const
{
    //XXX store 1. / (1. - rho_min)?
    const std::size_t Ai(particle_kinds_[i]);
    const real_type epsilon(std::get<2>(parameter_table_[Ai]));
    return (rho >= 1.)       ? 0. :
           (rho <= rho_min_) ? -c_scale_ * epsilon * c_linear_ :
           -c_scale_ * epsilon * (c_linear_ + 0.5 * (1. - c_linear_) *
           std::sin(constants<real_type>::pi * (1. - rho) / (1 - rho_min_)) *
           constants<real_type>::pi / (1. - rho_min_));
}


} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
