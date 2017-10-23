
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Optimization/missionSegmentSettings.h"

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_decision_variable_settings )

struct ObjectiveFunctionClass{

    ObjectiveFunctionClass( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > >
                            dynamicsSimulator, const double finalEpoch ) : dynamicsSimulator_( dynamicsSimulator ),
                                    finalEpoch_( finalEpoch ){ }

    ~ObjectiveFunctionClass( ){ }

    double userDefinedObjectiveFunction()
    {
        std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                dynamicsSimulator_->getEquationsOfMotionNumericalSolution().rbegin();

        double time1 = rit->first;
        Eigen::Vector6d state1 = rit->second.segment(0, 6);

        ++rit;

        double time2 = rit->first;
        Eigen::Vector6d state2 = rit->second.segment(0, 6);

        Eigen::Vector6d finalState = state2 + ( state1- state2)*(finalEpoch_ - time2)/(time1-time2);

        Eigen::Vector6d finalSphericalState =
                orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                    finalState);

        double radius = finalSphericalState[ orbital_elements::radius ];
        double flightPathAngle = finalSphericalState[ orbital_elements::flightPath ];
        double globalLongitude = finalSphericalState[ orbital_elements::longitude ];

        return (radius - 42000) + flightPathAngle*10000 + (globalLongitude - 1.41)*10000;

    }

    boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator_;
    double finalEpoch_;

};

BOOST_AUTO_TEST_CASE( testMisionSegmentSettings )
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

    boost::shared_ptr< SingleDecisionVariableSettings > decisionVarSettings1 = boost::make_shared<
            SingleDecisionVariableSettings >( simulation_time_decision_variable, 4500, 5000 );

    boost::shared_ptr< SingleDecisionVariableSettings > decisionVarSettings2 =
            boost::make_shared< SingleSphericalOrbitalElementDecisionVariableSettings >(
                    orbital_elements::speed, 9000, 10000, "SC1" );

    std::vector< boost::shared_ptr< SingleDecisionVariableSettings > >
            multipleDecisionVariableSettings;

    multipleDecisionVariableSettings.push_back( decisionVarSettings1 );
    multipleDecisionVariableSettings.push_back( decisionVarSettings2 );

    boost::shared_ptr< DecisionVariableSettings > decisionVariableSettings =
            boost::make_shared< DecisionVariableSettings >( multipleDecisionVariableSettings );

    boost::shared_ptr< ObjectiveFunctionClass > objectiveFunctionObject =
            boost::make_shared< ObjectiveFunctionClass >( dynamicsSimulator, 1000 );

    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings =
            boost::make_shared< UserDefinedObjectiveFunctionSettings >(
                boost::bind( &ObjectiveFunctionClass::userDefinedObjectiveFunction, objectiveFunctionObject ) );

    MissionSegmentSettings testBench1( dynamicsSimulator, decisionVariableSettings,
                                       objectiveFunctionSettings );

    MissionSegmentSettings testBench2( dynamicsSimulator );

    MissionSegmentSettings testBench3( dynamicsSimulator, decisionVarSettings1 );

    MissionSegmentSettings testBench4( dynamicsSimulator, objectiveFunctionSettings );

    dynamicsSimulator->integrateEquationsOfMotion( initialStates );

    BOOST_CHECK_EQUAL( testBench4.getFinalSimulationEpoch(), 1000.0 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
