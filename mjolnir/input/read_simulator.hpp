#ifndef MJOLNIR_READ_SIMULATOR
#define MJOLNIR_READ_SIMULATOR
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MDSimulator.hpp>
#include <mjolnir/util/make_unique.hpp>
#include "read_system.hpp"
#include "read_forcefield.hpp"
#include "read_integrator.hpp"
#include "read_observer.hpp"

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator(const toml::Table& data)
{
    const auto& simulator = data.at("simulator").cast<toml::value_t::Table>();
    const std::string type = toml::get<std::string>(simulator.at("simulator_type"));
    if(type == "Molecular Dynamics")
    {//XXX: separate this block from read_simulator() ?
        const std::string integration =
            toml::get<std::string>(simulator.at("time_integration_scheme"));
        const std::size_t tstep =
            toml::get<std::size_t>(simulator.at("total_step"));

        if(integration == "Newtonian")
        {
            return make_unique<
                MDSimulator<traitsT, VelocityVerletStepper<traitsT>>>(tstep,
                    read_system<traitsT>(data, 0), read_forcefield<traitsT>(data, 0),
                    read_velocity_verlet_stepper<traitsT>(data),
                    read_observer<traitsT>(data));
        }
        else if(integration == "Underdamped Langevin")
        {
            return make_unique<
                MDSimulator<traitsT, UnderdampedLangevinStepper<traitsT>>>(tstep,
                    read_system<traitsT>(data, 0), read_forcefield<traitsT>(data, 0),
                    read_underdamped_langevin_stepper<traitsT>(data),
                    read_observer<traitsT>(data));
        }
        else
        {
            throw std::runtime_error("invalid integration scheme: " + integration);
        }
    }
    else
    {
        throw std::runtime_error("invalid simulator type: " + type);
    }
}

}
#endif// MJOLNIR_READ_SIMULATOR
