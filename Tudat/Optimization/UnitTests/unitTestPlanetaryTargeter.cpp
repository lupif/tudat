
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Optimization/simplePlanetaryTargeter.h"

#include "pagmo/algorithm.hpp"

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_planetary_targeter )

BOOST_AUTO_TEST_CASE( testSimplePlanetaryTargeter )
{

    pagmo::random_device::set_seed( 123456 );

    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "naif0009.tls" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    double departureEpoch = spice_interface::convertDateStringToEphemerisTime( "2013-11-18T18:28:00" );

    double arrivalEpoch = spice_interface::convertDateStringToEphemerisTime( "2014-09-22T02:24:00" );

    optimization::UnperturbedNumericalPlanetaryTargeter planetTargeter( "Earth", "Mars", departureEpoch,
                                                                        arrivalEpoch, 3600.0 );


    std::string filePath_( __FILE__ );

    std::string outputFolder = filePath_.substr( 0, filePath_.length( ) -
            std::string( "unitTestPlanetaryTargeter.cpp" ).length( ) )  + "planetaryTargeterTestResults";

    input_output::writeDataMapToTextFile( planetTargeter.departurePlanetStateMap_,
                            "Earth.dat",
                            outputFolder,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    input_output::writeDataMapToTextFile( planetTargeter.spacecraftStateMap_,
                            "SC.dat",
                            outputFolder,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    input_output::writeDataMapToTextFile( planetTargeter.targetPlanetStateMap_,
                            "Mars.dat",
                            outputFolder,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

}


BOOST_AUTO_TEST_SUITE_END( )

}

}
