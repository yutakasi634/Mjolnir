#define BOOST_TEST_MODULE "test_directional_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/DirectionalContactInteraction.hpp>
#include <mjolnir/potential/local/CosinePotential.hpp>
#include <mjolnir/potential/local/GoContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/math/Matrix.hpp>

#include <tuple>

BOOST_AUTO_TEST_CASE(DirectionalContactInteraction_numerical_diff)
{
  mjolnir::LoggerManager::set_default_logger("test_directional_contact_interaction.log");
  using traits_type            = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
  using real_type              = traits_type::real_type;
  using coord_type             = traits_type::coordinate_type;
  using boundary_type          = traits_type::boundary_type;
  using system_type            = mjolnir::System<traits_type>;
  using angle_potential_type   = mjolnir::CosinePotential<real_type>;
  using contact_potential_type = mjolnir::GoContactPotential<real_type>;
  using directional_contact_type = mjolnir::DirectionalContactInteraction<
    traits_type, angle_potential_type, angle_potential_type, contact_potential_type>;

  const real_type k_angle(1e0);
  const real_type k_contact(1e0);
  const real_type angle1_native(0.0);
  const real_type angle2_native(mjolnir::math::constants<real_type>::pi() / 2.0);
  const real_type contact_native(2.0);
  angle_potential_type angle1_potential = angle_potential_type{k_angle, 1, angle1_native};
  angle_potential_type angle2_potential = angle_potential_type{k_angle, 2, angle2_native};
  contact_potential_type contact_potential = contact_potential_type{k_contact, contact_native};
  directional_contact_type interaction("none",
      {{std::make_tuple(std::array<std::size_t, 4>{0,1,2,3},
      angle1_potential, angle2_potential, contact_potential)
      }});

  system_type sys(4, boundary_type{});
  sys.mass(0) = 1.0;
  sys.mass(1) = 1.0;
  sys.mass(2) = 1.0;
  sys.mass(3) = 1.0;

  const int angle_step_num = 108;
  const std::size_t contact_step_num = 24;
  const real_type dtheta =
    mjolnir::math::constants<real_type>::pi() / (real_type(angle_step_num) / 2.0);
  const real_type contact_cutoff = contact_potential.cutoff();
  const real_type test_contact_range = 1.1 * contact_cutoff;
  const real_type dr = test_contact_range / real_type(contact_step_num);
  const real_type tol = 1e-4;
  const real_type dx = 1e-4;

  for(int i = - angle_step_num / 2; i < angle_step_num / 2; ++i)
    {
      for(int j = - angle_step_num / 2; j < angle_step_num / 2; ++j)
        {
          for(std::size_t k = 0; k < contact_step_num; ++k)
            {
              real_type theta1 = i * dtheta;
              real_type theta2 = j * dtheta;
              real_type r = k * dr;

              mjolnir::math::Matrix<real_type,3, 3> rotate1(
                  std::cos(theta1), -std::sin(theta1), 0.0,
                  std::sin(theta1), std::cos(theta1), 0.0,
                  0.0,              0.0,               1.0
                                                            );
              mjolnir::math::Matrix<real_type,3, 3> rotate2(
                  std::cos(theta2), -std::sin(theta2), 0.0,
                  std::sin(theta2), std::cos(theta2), 0.0,
                  0.0,              0.0,               1.0
                                                            );
              const coord_type x_unit_vec = {1.0, 0.0, 0.0};
              const coord_type particle1_vec = rotate1 * x_unit_vec;
              const coord_type particle2_vec = rotate2 * x_unit_vec;
              const coord_type r_vec = {r, 0, 0};
              const coord_type pos2(0e0, 0e0, 0e0);
              const coord_type pos1 = pos2 + particle1_vec;
              const coord_type pos4 = pos2 + r_vec;
              const coord_type pos3 = pos4 + particle2_vec;

              sys.at(0).position = pos1;
              sys.at(1).position = pos2;
              sys.at(2).position = pos3;
              sys.at(3).position = pos4;
              sys.at(0).velocity = coord_type(0.0, 0.0, 0.0);
              sys.at(1).velocity = coord_type(0.0, 0.0, 0.0);
              sys.at(2).velocity = coord_type(0.0, 0.0, 0.0);
              sys.at(3).velocity = coord_type(0.0, 0.0, 0.0);
              sys.at(0).force = coord_type(0.0, 0.0, 0.0);
              sys.at(1).force = coord_type(0.0, 0.0, 0.0);
              sys.at(2).force = coord_type(0.0, 0.0, 0.0);
              sys.at(3).force = coord_type(0.0, 0.0, 0.0);

              const auto init = sys;

              for(std::size_t idx=0; idx<4; ++idx)
                {
                  {
                    // reset positions
                    sys = init;

                    // calcu U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::X(sys.position(idx)) += dx;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::X(sys.position(idx)) += dx;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                  }
                  {
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Y(sys.position(idx)) += dx;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Y(sys.position(idx)) += dx;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    //central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                  }
                  {
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Z(sys.position(idx)) += dx;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Z(sys.position(idx)) += dx;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    //central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                  }
                }
            }
        }
    }
}
