#define BOOST_TEST_MODULE "test_read_uniform_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_local_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value env;
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
        )"_toml;

        const auto g = mjolnir::read_uniform_potential<real_type>(v, env);
        BOOST_TEST(g.k()  == real_type(3.14), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto env = u8R"(
            indices = [1, 2]
            k       = 3.14
        )"_toml;
        const auto v = u8R"(
            indices = "indices"
            k       = "k"
        )"_toml;
        const auto g = mjolnir::read_uniform_potential<real_type>(v, env);
        BOOST_TEST(g.k()  == real_type(3.14), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_uniform_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1,2], k = 3.14}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::UniformPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()  == real_type(3.14), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_uniform_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            env.indices = [1, 2]
            env.pi      = 3.14
            parameters = [
                {indices = "indices", k = "pi"}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::UniformPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()  == real_type(3.14), tolerance<real_type>());
    }
}
