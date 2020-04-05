#ifndef MJOLNIR_INPUT_READ_INPUT_FILE_HPP
#define MJOLNIR_INPUT_READ_INPUT_FILE_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/input/read_units.hpp>
#include <mjolnir/input/read_path.hpp>

#ifdef MJOLNIR_WITH_OPENMP
#include <mjolnir/omp/omp.hpp>
#endif

#include <memory>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
std::unique_ptr<SimulatorBase>
read_parallelism(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(simulator.as_table().count("parallelism") == 0)
    {
        MJOLNIR_LOG_NOTICE("execute on single core");
        return read_units<SimulatorTraits<realT, boundaryT>>(root, simulator);
    }

    const auto parallelism = toml::find(simulator, "parallelism");
    if(parallelism.is_string() && parallelism.as_string() == "sequencial")
    {
        MJOLNIR_LOG_NOTICE("execute on single core");
        return read_units<SimulatorTraits<realT, boundaryT>>(root, simulator);
    }
    else if(parallelism.is_string() &&
           (parallelism.as_string() == "openmp" ||
            parallelism.as_string() == "OpenMP"))
    {
#ifdef MJOLNIR_WITH_OPENMP
        MJOLNIR_LOG_NOTICE("execute on ", omp_get_max_threads() ," cores with openmp");
        return read_units<OpenMPSimulatorTraits<realT, boundaryT>>(root, simulator);
#else
        MJOLNIR_LOG_WARN("OpenMP flag is set, but OpenMP is not enabled when building.");
        MJOLNIR_LOG_WARN("Cannot use OpenMP, running with single core.");
        return read_units<SimulatorTraits<realT, boundaryT>>(root, simulator);
#endif
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_parallelism: invalid variable ",
            toml::find(simulator, "parallelism"), "here", {
            "- \"sequencial\": run with only 1 core (default)",
            "- \"openmp\"    : use openmp to parallelize."
            }));
    }
}

template<typename realT>
std::unique_ptr<SimulatorBase>
read_boundary(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto boundary = toml::find<std::string>(simulator, "boundary_type");
    if(boundary == "Unlimited")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Unlimited");
        return read_parallelism<realT, UnlimitedBoundary>(root, simulator);
    }
    else if(boundary == "Periodic")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is Periodic. "
                           "The shape is cuboid.");
        return read_parallelism<realT, CuboidalPeriodicBoundary>(root, simulator);
    }
    else if(boundary == "PeriodicCuboid")
    {
        MJOLNIR_LOG_NOTICE("Boundary Condition is PeriodicCuboid");
        return read_parallelism<realT, CuboidalPeriodicBoundary>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_boundary: invalid boundary",
            toml::find(simulator, "boundary_type"), "here", {
            "- \"Unlimited\"     : no boundary condition. infinite space",
            "- \"Periodic\"      : periodic boundary. Assuming cuboidal shape.",
            "- \"PeriodicCuboid\": periodic boundary with cuboidal shape"
            }));
    }
}

inline std::unique_ptr<SimulatorBase>
read_precision(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto prec = toml::find<std::string>(simulator, "precision");
    if(prec == "double")
    {
        MJOLNIR_LOG_NOTICE("precision is double");
        return read_boundary<double>(root, simulator);
    }
    else if(prec == "float")
    {
        MJOLNIR_LOG_NOTICE("precision is float");
        return read_boundary<float>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_precision: invalid precision",
            toml::find(simulator, "precition"), "here", {
            "expected value is one of the following.",
            "- \"double\": 64 bit floating-point",
            "- \"float\" : 32 bit floating-point"
            }));
    }
}

void try_merge_val(toml::value& first_val, toml::value& second_val,
                   const std::vector<toml::key>& subject_branch_keys,
                   const toml::string& include_file)
{
    toml::value old_first_val = first_val;
    if(!(old_first_val.is_table() && second_val.is_table()))
    {
        std::string branch_keys_str = "";
        for(const auto& key : subject_branch_keys)
        {
            branch_keys_str += "." + std::string(key);
        }
        branch_keys_str.erase(0,1);

        throw_exception<std::runtime_error>("[error] "
        "mjolnir::expand_include_file reading `", std::string(include_file), "`: ",
        branch_keys_str, " definition is duplicated");
    }
    else
    {
        for(auto& key_val_in_second : second_val.as_table())
        {
            auto& second_table_key = key_val_in_second.first;
            auto& second_table_val = key_val_in_second.second;
            if(old_first_val.contains(second_table_key))
            {
                auto new_subject_branch_keys = subject_branch_keys;
                new_subject_branch_keys.push_back(second_table_key);
                try_merge_val(toml::find(first_val, second_table_key),
                              second_table_val, new_subject_branch_keys, include_file);
            }
            else
            {
                first_val.as_table().insert(key_val_in_second);
            }
        }
    }
    return;
}

void expand_include_file(toml::value& root, const std::vector<toml::key>& subject_branch_keys)
{
    while(root.contains("include_file"))
    {
        const auto include_file = toml::find(root, "include_file");
        root.as_table().erase("include_file");
        if(include_file.is_array())
        {
            for(auto include_file_val : include_file.as_array())
            {
                std::string include_file_name = include_file_val.as_string();
                std::cerr << "-- expanding toml file `" << include_file_name << "` ... ";
                toml::value val_in_file = toml::parse(include_file_name);
                try_merge_val(root, val_in_file, subject_branch_keys, include_file_name);
                std::cerr << " successfully expanded." << std::endl;
            }
        }
        else
        {
            std::string include_file_name = include_file.as_string();
            std::cerr << "-- expanding toml file `" << include_file.as_string() << "` ... ";
            toml::value val_in_file = toml::parse(include_file_name);
            try_merge_val(root, val_in_file, subject_branch_keys, include_file_name);
            std::cerr << " successfully expanded." << std::endl;
        }
    }

    for(auto& key_value : root.as_table())
    {
        auto new_branch_keys = subject_branch_keys;
        new_branch_keys.push_back(key_value.first);
        if(key_value.second.is_table())
        {
            expand_include_file(key_value.second, new_branch_keys);
        }
        else if(key_value.second.is_array())
        {
            for(auto& value : key_value.second.as_array())
            {
                if(value.is_table())
                {
                    expand_include_file(value, new_branch_keys);
                }
            }
        }
    }
    return ;
}

inline std::unique_ptr<SimulatorBase>
read_input_file(const std::string& filename)
{
    // here, logger name is not given yet. output status directory on console.
    std::cerr << "-- reading and parsing toml file `" << filename << "` ... ";
    auto root = toml::parse(filename);
    std::cerr << " successfully parsed." << std::endl;

    // expand include_file key in toml value tree.
    expand_include_file(root, {toml::key("root")});

    // initializing logger by using output_path and output_prefix ...
    const auto& output   = toml::find(root, "files", "output");
    const auto  out_path = read_output_path(root);

    // XXX:  Here, this code assumes POSIX. it does not support windows.
    // TODO: Consider using Boost.filesystem to manage path and files
    //       in more elegant and powerful way? After switching C++17,
    //       we can re-write that with <filesystem>.
    const auto logger_name = out_path +
        toml::find<std::string>(output, "prefix") + ".log";

    MJOLNIR_SET_DEFAULT_LOGGER(logger_name);
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_NOTICE("the log file is `", logger_name, '`');

    // Check top-level toml-values. Since it uses logger to warn,
    // we need to call it after `MJOLNIR_SET_DEFAULT_LOGGER(logger_name)`.
    check_keys_available(root, {"files"_s, "units"_s, "simulator"_s,
                                "systems"_s, "forcefields"_s, "include_file"_s});

    // the most of important flags are defined in [simulator], like
    // `precision = "float"`, `boundary_type = "Unlimited"`.
    // Those values should be read before others.
    // Thus first read [simulator] here and pass it to the latter functions.

    const auto& simulator = toml::find(root, "simulator");
    return read_precision(root,  simulator);
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<SimulatorBase> read_parallelism<double, UnlimitedBoundary       >(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_parallelism<float , UnlimitedBoundary       >(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_parallelism<double, CuboidalPeriodicBoundary>(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_parallelism<float , CuboidalPeriodicBoundary>(const toml::value& root, const toml::value& simulator);

extern template std::unique_ptr<SimulatorBase> read_boundary<double>(const toml::value& root, const toml::value& simulator);
extern template std::unique_ptr<SimulatorBase> read_boundary<float >(const toml::value& root, const toml::value& simulator);
#endif

}// mjolnir
#endif// MJOLNIR_READ_INPUT_FILE
