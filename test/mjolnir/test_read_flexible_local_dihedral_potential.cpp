#define BOOST_TEST_MODULE "test_read_flexible_local_dihedral_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_potential.hpp>

BOOST_AUTO_TEST_CASE(read_flexible_local_dihedral_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_dihedral.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {"coef",    toml::value{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}},
        };
        const auto g = mjolnir::read_flexible_local_dihedral_potential<real_type>(v);
        BOOST_TEST(g.k()       == 3.14, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[0] ==  1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[1] ==  2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[2] ==  3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[3] ==  4.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[4] ==  5.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[5] ==  6.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[6] ==  7.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_flexible_local_dihedral_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_dihedral.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {"coef",    toml::value{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}},
        };
        const auto g = mjolnir::read_flexible_local_dihedral_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[0] ==  1.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[1] ==  2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[2] ==  3.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[3] ==  4.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[4] ==  5.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[5] ==  6.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[6] ==  7.0f, boost::test_tools::tolerance(tol));
    }
}
