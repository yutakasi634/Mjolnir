#ifndef MJOLNIR_FORCEFIELD_3SPN2_BOND_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_BOND_POTENTIAL_HPP
#include <limits>

namespace mjolnir
{

template<typename T> class System;

template<typename realT>
class ThreeSPN2BondPotential
{
  public:
    using real_type = realT;

  public:
    ThreeSPN2BondPotential(const real_type k, const real_type v0)
        : k_(k), v0_(v0)
    {}
    ~ThreeSPN2BondPotential() = default;

    // k * (v - v0)^2 + 100k * (v - v0)^4
    real_type potential(const real_type v) const noexcept
    {
        const real_type dv  = v - v0_;
        const real_type dv2 = dv * dv;
        return this->k_ * (dv2 + real_type(100) * dv2 * dv2);
    }

    // 2k * (v - v0) + 400k * (v - v0)^3
    real_type derivative(const real_type v) const noexcept
    {
        const real_type dv = v - v0_;
        return 2 * this->k_ * (dv + 200 * dv * dv * dv);
    }

    template<typename T>
    void initialize(const System<T>&) const noexcept {return;}

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "ThreeSPN2BondPotential";}

    real_type k()  const noexcept {return k_;}
    real_type v0() const noexcept {return v0_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:

    real_type k_;
    real_type v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ThreeSPN2BondPotential<double>;
extern template class ThreeSPN2BondPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_HARMONIC_POTENTIAL */
