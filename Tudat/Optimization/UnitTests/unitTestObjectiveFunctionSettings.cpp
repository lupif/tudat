
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Optimization/objectiveFunction.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"


struct UserDefinedStructure
{
    UserDefinedStructure( boost::shared_ptr<
                          tudat::propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator )
    : dynamicsSimulator_( dynamicsSimulator ){ }

    ~UserDefinedStructure( ){ }

    double theUserDefinedFunction( )
    {
        std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                dynamicsSimulator_->getEquationsOfMotionNumericalSolution().rbegin();
        return( tudat::orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                    rit->second.segment( 0, 6 ) )[tudat::orbital_element_conversions::speedIndex] );
    }

    boost::shared_ptr< tudat::propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator_;
};

namespace tudat
{

namespace unit_tests
{

using namespace optimization;



BOOST_AUTO_TEST_SUITE( test_decision_variable_settings )

BOOST_AUTO_TEST_CASE( testObjectiveFunctionConstructors )
{
    using namespace optimization;
    using namespace propagators;

    boost::shared_ptr< SingleDependentVariableSaveSettings > depVarSaveSettings =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                propagators::relative_distance_dependent_variable, "SC1", "Earth" );

    ObjectiveFunctionFromFinalDependentVariableSettings
            objectiveFunctionSettings( depVarSaveSettings, 7000000.0 );

    fflush( stdout );
}

BOOST_AUTO_TEST_CASE( testObjectiveFunctionMethods )
{

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
    bodiesToCreate.push_back( "Moon" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, 0.0, 100300.0 );

    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    bodySettings[ "Moon" ]->rotationModelSettings->resetOriginalFrame( "J2000" );


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
                    relative_distance_dependent_variable, "SC1", "SC2") );
    depVarSaves.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                    relative_speed_dependent_variable, "SC1", "SC2") );

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


    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationSettings,
                    bodiesToPropagate, initialStates, 1000, cowell,
                            depVarsaveSettings);


    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared<
            IntegratorSettings< > >(
                rungeKutta4, 0.0, 10 );

    boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator = boost::make_shared<
            SingleArcDynamicsSimulator< > >(
                bodyMap, integratorSettings, propagatorSettings, false );

    boost::shared_ptr< UserDefinedStructure > usrDefinedStructure
            = boost::make_shared< UserDefinedStructure >( dynamicsSimulator );

    dynamicsSimulator->integrateEquationsOfMotion( initialStates );

    boost::function< double () > userDefinedFuction = boost::bind(
                &UserDefinedStructure::theUserDefinedFunction, usrDefinedStructure );

    UserDefinedObjectiveFunctionSettings userDefinedObjectiveFunction( userDefinedFuction );

    std::map< double, Eigen::VectorXd >::reverse_iterator rit =
            dynamicsSimulator->getEquationsOfMotionNumericalSolution().rbegin();

    double finalSpeed = convertCartesianToSphericalOrbitalState(
                rit->second.segment( 0, 6 ) )[speedIndex];

    BOOST_CHECK_EQUAL( userDefinedObjectiveFunction.userDefinedFunction_(), finalSpeed );

    FinalSphericalOrbitalElementObjectiveFunctionSettings testBench2( radius, "SC1", 0.0 );

    ObjectiveFunctionFromFinalDependentVariableSettings testBench3( boost::make_shared< SingleDependentVariableSaveSettings >(
        relative_speed_dependent_variable, "SC1", "SC2" ), 0.0 );

    BOOST_CHECK_EQUAL( testBench3.getPositionInDependentVariablesSaveVector( dynamicsSimulator ), 1 );

    ObjectiveFunctionFromMinOrMaxDependentVariableSettings testBench4( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                                           relative_speed_dependent_variable, "SC1", "SC2" ), 0.0 );

    ObjectiveFunctionFromMinOrMaxDependentVariableSettings testBench5( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                                           relative_distance_dependent_variable, "SC1", "SC2" ), 0.0,
                                                                       0, false );

    std::map< double, Eigen::VectorXd > dependentVariableMap = dynamicsSimulator->getDependentVariableHistory();
    std::map< double, Eigen::VectorXd >::iterator it = dependentVariableMap.begin();

    double maxDistance = it->second[0];
    double minSpeed = it->second[1];

    while( it != dependentVariableMap.end() )
    {

        if( maxDistance < it->second[0] )
            maxDistance = it->second[0];
        if( minSpeed > it->second[1] )
            minSpeed = it->second[1];
        ++it;
    }

    BOOST_CHECK_EQUAL( testBench4.getMinOrMaxSaveVariable( dynamicsSimulator ), minSpeed );

    BOOST_CHECK_EQUAL( testBench5.getMinOrMaxSaveVariable( dynamicsSimulator ), maxDistance );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
