#include <mjolnir/input/read_global_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ExcludedVolumePotential<double> read_excluded_volume_potential(const toml::value& global);
template ExcludedVolumePotential<float > read_excluded_volume_potential(const toml::value& global);

template HardCoreExcludedVolumePotential<double> read_hard_core_excluded_volume_potential(const toml::value& global);
template HardCoreExcludedVolumePotential<float > read_hard_core_excluded_volume_potential(const toml::value& global);

template LennardJonesPotential<double> read_lennard_jones_potential(const toml::value& global);
template LennardJonesPotential<float > read_lennard_jones_potential(const toml::value& global);

template UniformLennardJonesPotential<double> read_uniform_lennard_jones_potential(const toml::value& global);
template UniformLennardJonesPotential<float > read_uniform_lennard_jones_potential(const toml::value& global);

template DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_debye_huckel_potential(const toml::value& global);

template ThreeSPN2ExcludedVolumePotential<double> read_3spn2_excluded_volume_potential(const toml::value& global);
template ThreeSPN2ExcludedVolumePotential<float > read_3spn2_excluded_volume_potential(const toml::value& global);
}
