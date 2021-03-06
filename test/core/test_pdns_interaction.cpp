#define BOOST_TEST_MODULE "test_PDNS_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificInteraction.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>

BOOST_AUTO_TEST_CASE(PDNS_Interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using real_type              = double;
    using potential_type         = mjolnir::ProteinDNANonSpecificPotential<traits_type>;
    using contact_parameter_type = typename potential_type::contact_parameter_type;
    using dna_index_type         = typename potential_type::dna_index_type;
    using interaction_type       = mjolnir::ProteinDNANonSpecificInteraction<traits_type>;
    using ignore_molecule_type   = typename potential_type::ignore_molecule_type;
    using ignore_group_type      = typename potential_type::ignore_group_type;
    using partition_type         = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    //  0 o              |
    //     \             |
    //    1 o === o 3    |
    //     /       \     |
    //  2 o         o 4  |

    constexpr auto pi  = mjolnir::math::constants<real_type>::pi();

    const real_type r0     = 5.0;
    const real_type theta0 = pi * 0.5;
    const real_type phi0   = pi * 2.0 / 3.0;

    const real_type sigma  = 1.0;
    const real_type delta  = pi / 18.0;

    const contact_parameter_type p_pro{1, 0, 2, -1.2, r0, theta0, phi0, r0+5.0, (r0+5.0)*(r0+5.0)};
    const dna_index_type p_dna{3, 4};

    potential_type potential(
        sigma, delta, 5.0, {p_pro}, {p_dna}, {/*exclusion*/},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    interaction_type interaction(potential_type(potential),
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(5, boundary_type{});
    topology_type topol(5);

    topol.add_connection(0, 1, "bond");
    topol.add_connection(1, 2, "bond");
    topol.add_connection(3, 4, "bond");
    topol.construct_molecules();

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)  = 1.0;
        sys.rmass(i) = 1.0;
        sys.position(i) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(i) = coord_type(0.0, 0.0, 0.0);
        sys.force(i)    = coord_type(0.0, 0.0, 0.0);
        sys.name(i)  = "DNA";
        sys.group(i) = "DNA";
    }
    potential  .initialize(sys, topol);
    interaction.initialize(sys, topol);

    for(const auto r     : {r0 - 0.2, r0 + 0.5})
    {
    for(const auto theta : {theta0 + delta * 0.2,
                            theta0 + delta * 0.7,
                            theta0 + delta * 1.2,
                            theta0 - delta * 1.2,
                            theta0 - delta * 0.7,
                            theta0 - delta * 0.2})
    {
    for(const auto phi   : {phi0 + delta * 0.2,
                            phi0 + delta * 0.7,
                            phi0 + delta * 1.2,
                            phi0 - delta * 1.2,
                            phi0 - delta * 0.7,
                            phi0 - delta * 0.2})
    {
        BOOST_TEST_MESSAGE("=======================================================");
        BOOST_TEST_MESSAGE("r = " << r << ", theta = " << theta << ", phi = " << phi);
        BOOST_TEST_MESSAGE("theta0 = " << theta0 << ", delta = " << delta);
        //  2 o              |
        //     \             |
        //  ^ 1 o === o 3    |
        //     /       \     |
        //  0 o         o 4  |


        sys.position(1) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
        sys.position(0) = mjolnir::math::make_coordinate<coord_type>(
                -std::cos(theta) - 1.0, -std::sin(theta), 0.0);
        sys.position(2) = mjolnir::math::make_coordinate<coord_type>(
                 std::cos(theta) - 1.0,  std::sin(theta), 0.0);

        sys.position(3) = mjolnir::math::make_coordinate<coord_type>(
                r, 0.0, 0.0);
        sys.position(4) = mjolnir::math::make_coordinate<coord_type>(
                r + std::cos(pi-phi), -std::sin(pi-phi), 0.0);

        // rotate in random direction
        {
            using matrix33_type = typename traits_type::matrix33_type;

            const auto rot_x = uni(mt) * pi;
            const auto rot_y = uni(mt) * pi;
            const auto rot_z = uni(mt) * pi;

            matrix33_type rotm_x(1.0,             0.0,              0.0,
                                 0.0, std::cos(rot_x), -std::sin(rot_x),
                                 0.0, std::sin(rot_x),  std::cos(rot_x));
            matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                              0.0, 1.0,              0.0,
                                 -std::sin(rot_y), 0.0,  std::cos(rot_y));
            matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                                 std::sin(rot_z),  std::cos(rot_z), 0.0,
                                             0.0,              0.0, 1.0);

            const matrix33_type rotm = rotm_x * rotm_y * rotm_z;
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                sys.position(idx) = rotm * sys.position(idx);
            }
        }
        // check configuration is okay
        {
            // XXX here it does not use boundary condition.
            const auto v13 = sys.position(3) - sys.position(1);
            BOOST_TEST_REQUIRE(mjolnir::math::length(v13) == r,
                               boost::test_tools::tolerance(1e-3));

            const auto v02 = sys.position(2) - sys.position(0);
            const auto dot_1302 = mjolnir::math::dot_product(v02, v13);
            const auto cos_1302 = mjolnir::math::clamp<real_type>(dot_1302 *
                mjolnir::math::rlength(v13) * mjolnir::math::rlength(v02), -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos_1302) == theta,
                               boost::test_tools::tolerance(1e-3));

            const auto v43 = sys.position(3) - sys.position(4);
            const auto dot_134 = mjolnir::math::dot_product(v13, v43);
            const auto cos_134 = mjolnir::math::clamp<real_type>(dot_134 *
                mjolnir::math::rlength(v13) * mjolnir::math::rlength(v43), -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos_134) == phi,
                               boost::test_tools::tolerance(1e-3));
        }

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
        }
        const system_type init = sys;
        interaction.initialize(init, topol);

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::X(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::X(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::X(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
                BOOST_TEST_MESSAGE("force.x = " << mjolnir::math::X(sys.force(idx)) <<
                                   ", -dE/dr = " << -dE/dr);
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::Y(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Y(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Y(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
                BOOST_TEST_MESSAGE("force.y = " << mjolnir::math::Y(sys.force(idx)) <<
                                   ", -dE/dr = " << -dE/dr);
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::Z(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Z(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Z(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
                BOOST_TEST_MESSAGE("force.z = " << mjolnir::math::Z(sys.force(idx)) <<
                                   ", -dE/dr = " << -dE/dr);
            }
        }
    } // theta2
    } // theta1
    } // rbp
}

BOOST_AUTO_TEST_CASE(PDNS_Interaction_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using real_type              = double;
    using potential_type         = mjolnir::ProteinDNANonSpecificPotential<traits_type>;
    using contact_parameter_type = typename potential_type::contact_parameter_type;
    using dna_index_type         = typename potential_type::dna_index_type;
    using interaction_type       = mjolnir::ProteinDNANonSpecificInteraction<traits_type>;
    using ignore_molecule_type   = typename potential_type::ignore_molecule_type;
    using ignore_group_type      = typename potential_type::ignore_group_type;
    using partition_type         = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    //  0 o              |
    //     \             |
    //    1 o === o 3    |
    //     /       \     |
    //  2 o         o 4  |

    constexpr auto pi  = mjolnir::math::constants<real_type>::pi();

    const real_type r0     = 5.0;
    const real_type theta0 = pi * 0.5;
    const real_type phi0   = pi * 2.0 / 3.0;

    const real_type sigma  = 1.0;
    const real_type delta  = pi / 18.0;

    const contact_parameter_type p_pro{1, 0, 2, -1.2, r0, theta0, phi0, r0+5.0, (r0+5.0)*(r0+5.0)};
    const dna_index_type p_dna{3, 4};

    potential_type potential(
        sigma, delta, 5.0, {p_pro}, {p_dna}, {/*exclusion*/},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    interaction_type interaction(potential_type(potential),
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(5, boundary_type{});
    topology_type topol(5);

    topol.add_connection(0, 1, "bond");
    topol.add_connection(1, 2, "bond");
    topol.add_connection(3, 4, "bond");
    topol.construct_molecules();

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)  = 1.0;
        sys.rmass(i) = 1.0;
        sys.position(i) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(i) = coord_type(0.0, 0.0, 0.0);
        sys.force(i)    = coord_type(0.0, 0.0, 0.0);
        sys.name(i)  = "DNA";
        sys.group(i) = "DNA";
    }
    potential  .initialize(sys, topol);
    interaction.initialize(sys, topol);

    for(const auto r     : {r0 - 0.2, r0 + 0.5})
    {
    for(const auto theta : {theta0 + delta * 0.2,
                            theta0 + delta * 0.7,
                            theta0 + delta * 1.2,
                            theta0 - delta * 1.2,
                            theta0 - delta * 0.7,
                            theta0 - delta * 0.2})
    {
    for(const auto phi   : {phi0 + delta * 0.2,
                            phi0 + delta * 0.7,
                            phi0 + delta * 1.2,
                            phi0 - delta * 1.2,
                            phi0 - delta * 0.7,
                            phi0 - delta * 0.2})
    {
        BOOST_TEST_MESSAGE("=======================================================");
        BOOST_TEST_MESSAGE("r = " << r << ", theta = " << theta << ", phi = " << phi);
        BOOST_TEST_MESSAGE("theta0 = " << theta0 << ", delta = " << delta);
        //  2 o              |
        //     \             |
        //  ^ 1 o === o 3    |
        //     /       \     |
        //  0 o         o 4  |


        sys.position(1) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
        sys.position(0) = mjolnir::math::make_coordinate<coord_type>(
                -std::cos(theta) - 1.0, -std::sin(theta), 0.0);
        sys.position(2) = mjolnir::math::make_coordinate<coord_type>(
                 std::cos(theta) - 1.0,  std::sin(theta), 0.0);

        sys.position(3) = mjolnir::math::make_coordinate<coord_type>(
                r, 0.0, 0.0);
        sys.position(4) = mjolnir::math::make_coordinate<coord_type>(
                r + std::cos(pi-phi), -std::sin(pi-phi), 0.0);

        // rotate in random direction
        {
            using matrix33_type = typename traits_type::matrix33_type;

            const auto rot_x = uni(mt) * pi;
            const auto rot_y = uni(mt) * pi;
            const auto rot_z = uni(mt) * pi;

            matrix33_type rotm_x(1.0,             0.0,              0.0,
                                 0.0, std::cos(rot_x), -std::sin(rot_x),
                                 0.0, std::sin(rot_x),  std::cos(rot_x));
            matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                              0.0, 1.0,              0.0,
                                 -std::sin(rot_y), 0.0,  std::cos(rot_y));
            matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                                 std::sin(rot_z),  std::cos(rot_z), 0.0,
                                             0.0,              0.0, 1.0);

            const matrix33_type rotm = rotm_x * rotm_y * rotm_z;
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                sys.position(idx) = rotm * sys.position(idx);
            }
        }
        // check configuration is okay
        {
            // XXX here it does not use boundary condition.
            const auto v13 = sys.position(3) - sys.position(1);
            BOOST_TEST_REQUIRE(mjolnir::math::length(v13) == r,
                               boost::test_tools::tolerance(1e-3));

            const auto v02 = sys.position(2) - sys.position(0);
            const auto dot_1302 = mjolnir::math::dot_product(v02, v13);
            const auto cos_1302 = mjolnir::math::clamp<real_type>(dot_1302 *
                mjolnir::math::rlength(v13) * mjolnir::math::rlength(v02), -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos_1302) == theta,
                               boost::test_tools::tolerance(1e-3));

            const auto v43 = sys.position(3) - sys.position(4);
            const auto dot_134 = mjolnir::math::dot_product(v13, v43);
            const auto cos_134 = mjolnir::math::clamp<real_type>(dot_134 *
                mjolnir::math::rlength(v13) * mjolnir::math::rlength(v43), -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos_134) == phi,
                               boost::test_tools::tolerance(1e-3));
        }

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
        }
        system_type ref_sys = sys;
        interaction.initialize(sys, topol);

        constexpr real_type tol = 1e-4;

        const auto energy = interaction.calc_force_and_energy(sys);
        const auto ref_energy = interaction.calc_energy(ref_sys);
        interaction.calc_force(ref_sys);
        BOOST_TEST(ref_energy == energy, boost::test_tools::tolerance(tol));

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            BOOST_TEST(mjolnir::math::X(sys.force(idx)) == mjolnir::math::X(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(sys.force(idx)) == mjolnir::math::Y(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(sys.force(idx)) == mjolnir::math::Z(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
        }

    } // theta2
    } // theta1
    } // rbp
}
