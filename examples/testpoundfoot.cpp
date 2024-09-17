#include <boost/units/base_units/us/foot.hpp>
#include <boost/units/base_units/us/pound_force.hpp>
#include <boost/units/io.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/prefixes.hpp>
#include <boost/units/systems/si/pressure.hpp>
#include <boost/units/systems/si/torque.hpp>
#include <iostream>

namespace boost {
namespace units {
namespace us {

typedef make_system<foot_base_unit, pound_force_base_unit>::type system;
typedef unit<torque_dimension, system> torque;

BOOST_UNITS_STATIC_CONSTANT(pound_feet, torque);

}   // namespace us
}   // namespace units
}   // namespace boost

using namespace boost::units;
using namespace boost::units::si;

int
main() {
    quantity<us::torque> colonial_measurement(1.0 * us::pound_feet);
    std::cerr << quantity<si::torque>(colonial_measurement) << std::endl;

    quantity<pressure> p(101.1 * kilo * pascals);
    double dblP = p / pascals;   // double value in Pascals
    std::cout << "With units: " << p << std::endl;
    std::cout << "Without units: " << dblP << std::endl;
    return 0;
}