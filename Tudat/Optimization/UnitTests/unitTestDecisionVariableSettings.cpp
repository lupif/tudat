
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Optimization/decisionVariableSettings.h"

namespace tudat
{
namespace unit_tests
{

using namespace optimization;

BOOST_AUTO_TEST_SUITE( test_decision_variable_settings )

BOOST_AUTO_TEST_CASE( testSingleDecisionVariableSettings )
{

    Eigen::VectorXd lowerBoundary;
    lowerBoundary.resize(3);
    lowerBoundary[0] = 20;
    lowerBoundary[1] = 20;
    lowerBoundary[2] = 20;

    Eigen::VectorXd upperBoundary;
    upperBoundary.resize(3);
    upperBoundary[0] = 25;
    upperBoundary[1] = 25;
    upperBoundary[2] = 25;

    SingleDecisionVariableSettings newTestBench(
                    initialVelocitytateCartesianComponents, lowerBoundary, upperBoundary );

}

BOOST_AUTO_TEST_CASE( testDecisionVariableFromTerminationSettings )
{

    using namespace tudat;
    using namespace optimization;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace orbital_element_conversions;
    using namespace orbital_elements;
    using namespace numerical_integrators;


    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, 0.0, 100300.0 );

    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "SC1" ] = boost::make_shared< Body >();
    bodyMap[ "SC2" ] = boost::make_shared< Body >();

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "SC1" );
    bodiesToPropagate.push_back( "SC2" );

    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Earth" );
    centralBodies.push_back( "Earth" );

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >
            accelerationSettingsSC1, accelerationSettingsSC2;
    accelerationSettingsSC1[ "Earth" ].push_back(
                boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationSettingsSC2[ "Earth" ].push_back(
                boost::make_shared< AccelerationSettings >( central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "SC1" ] = accelerationSettingsSC1;
    accelerationMap[ "SC2" ] = accelerationSettingsSC2;


    AccelerationMap accelerationSettings = simulation_setup::createAccelerationModelsMap( bodyMap, accelerationMap,
        bodiesToPropagate, centralBodies );

    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > depVarSaves;

    depVarSaves.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    relative_distance_dependent_variable, "SC1", "Earth") );
    depVarSaves.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                    relative_speed_dependent_variable, "SC1", "Earth") );

    Eigen::VectorXd initialStates;

    initialStates.resize(12);

    Eigen::Vector6d initialKeplerElements = Eigen::Vector6d::Zero();
    initialKeplerElements[ semiMajorAxisIndex ] = 20000000;
    initialKeplerElements[ eccentricityIndex ] = 0.2;
    initialKeplerElements[ inclinationIndex ] = 20*mathematical_constants::PI/180;

    initialStates.segment(0, 6) = convertKeplerianToCartesianElements( initialKeplerElements,
            bodyMap[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() );

    initialKeplerElements[ trueAnomalyIndex ] = mathematical_constants::PI/3;

    initialStates.segment(6, 6) = convertKeplerianToCartesianElements( initialKeplerElements,
            bodyMap[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() );


    boost::shared_ptr< DependentVariableSaveSettings > depVarsaveSettings =
            boost::make_shared< DependentVariableSaveSettings >(depVarSaves, false);

    std::vector< boost::shared_ptr< PropagationTerminationSettings > > vecOfTerminationSettings;

    vecOfTerminationSettings.push_back( boost::make_shared< PropagationTimeTerminationSettings >(
                                            1000 ) );
    vecOfTerminationSettings.push_back( boost::make_shared<
                                        PropagationDependentVariableTerminationSettings >(
                                            depVarSaves[0], 1000000, false ) );

    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationHybridTerminationSettings >( vecOfTerminationSettings );

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationSettings,
                    bodiesToPropagate, initialStates, terminationSettings, cowell,
                            depVarsaveSettings);

    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared<
            IntegratorSettings< > >(
                rungeKutta4, 0.0, 10 );


    boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator = boost::make_shared<
            SingleArcDynamicsSimulator< > >(
                bodyMap, integratorSettings, propagatorSettings, false );

    double lowerBoundary = 800;
    double upperBoundary = 900;

    SingleDecisionVariableFromTerminationSettings testBench1( dynamicsSimulator, lowerBoundary,
                                                              upperBoundary );

    BOOST_CHECK_EQUAL( *(testBench1.memoryPositionOfVariable_),
                       1000 );

    BOOST_CHECK_EQUAL(testBench1.decisionVariable_,
                      fromTerminationSettings );

    *(testBench1.memoryPositionOfVariable_) = 2000;

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsRetrieve =
            boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< > >(
                dynamicsSimulator->getPropagatorSettings() );
    boost::shared_ptr< PropagationHybridTerminationSettings > terminationSettingsRetrieve =
            boost::dynamic_pointer_cast< PropagationHybridTerminationSettings >(
                propagatorSettingsRetrieve->getTerminationSettings());

    BOOST_CHECK_EQUAL(boost::dynamic_pointer_cast< PropagationTimeTerminationSettings >(
            terminationSettingsRetrieve->terminationSettings_[0])->terminationTime_, 2000 );

    SingleKeplerElementDecisionVariableSettings testBench2( inclination, "Earth", -5, 5, "SC2" );

    //modifyKeplerElement( inclination, dynamicsSimulator, "SC2", "Earth", mathematical_constants::PI/4 );

    Eigen::Vector6d modifiedCartState = propagatorSettingsRetrieve->getInitialStates().segment(6,6);

    //Eigen::Vector6d modifiedKeplerElements = convertCartesianToKeplerianElements( modifiedCartState,
    //        bodyMap["Earth"]->getGravityFieldModel()->getGravitationalParameter() );

    //BOOST_CHECK_CLOSE(mathematical_constants::PI/4, modifiedKeplerElements[inclinationIndex], 1e-3);

    SingleSphericalOrbitalElementDecisionVariableSettings testBench3( speed, 8000, 8500, "SC1" );

    //modifySphericalComponent( speed, dynamicsSimulator, "SC1", 8250.0 );

    //modifiedCartState = propagatorSettingsRetrieve->getInitialStates().segment(0,6);

    //Eigen::Vector6d modifiedSphericalElements = convertCartesianToSphericalOrbitalState( modifiedCartState );

    //BOOST_CHECK_CLOSE( modifiedSphericalElements[ speedIndex ], 8250.0, 1e-3);

    //SingleCartesianComponentDecisionVariableSettings testBench4( zCartesianVelocity, -200, 200, "SC2" );

    //modifyCartesianComponent( zCartesianVelocity, dynamicsSimulator, "SC2", 100.0 );

    //BOOST_CHECK_EQUAL( propagatorSettingsRetrieve->getInitialStates().segment(6,6)[5], 100.0 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
