#ifndef MJOLNIR_INPUT_READ_PATH_HPP
#define MJOLNIR_INPUT_READ_PATH_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this function may be called from other read_* functions.
template<typename C, template<typename...> class T, template<typename...> class A>
std::string read_input_path(const toml::basic_value<C, T, A>& root)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto& files = toml::find(root, "files");

    std::string input_path("./");
    if(files.contains("input"))
    {
        const auto& input = toml::find(files, "input");
        if(input.contains("path"))
        {
            input_path = toml::find<std::string>(input, "path");
            if(input_path.back() != '/') {input_path += '/';}
        }
    }
    return input_path;
}

template<typename C, template<typename...> class T, template<typename...> class A>
std::string read_output_path(const toml::basic_value<C, T, A>& root)
{
    // This function does not output log because it might be called when
    // the logger has not been initialized yet.
    // Logging should be done after calling this function.

    const auto& files = toml::find(root, "files");

    std::string output_path("./");
    if(files.contains("output"))
    {
        const auto& output = toml::find(files, "output");
        if(output.contains("path"))
        {
            output_path = toml::find<std::string>(output, "path");
            if(output_path.back() != '/') {output_path += '/';}
        }
    }
    return output_path;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::string read_input_path <toml::discard_comments, std::unordered_map, std::vector>(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>&);
extern template std::string read_output_path<toml::discard_comments, std::unordered_map, std::vector>(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>&);
#endif

} // mjolnir
#endif// MJOLNIR_READ_PATH_HPP
