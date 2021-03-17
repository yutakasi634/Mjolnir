#include <mjolnir/omp/PeriodicGridStoichiometricInteractionCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GlobalStoichiometricInteractionPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GlobalStoichiometricInteractionPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
