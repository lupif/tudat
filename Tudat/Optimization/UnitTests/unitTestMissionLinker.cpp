
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Optimization/missionLinker.h"
#include "pagmo/island.hpp"


namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_mission_linker )



struct ObjectiveFunction{

    ObjectiveFunction(  boost::shared_ptr< optimization::MissionSegmentSettings > missionSegment,
                        const double desiredLongitude, const double desiredRadius)
        : missionSegment_(missionSegment), desiredLongitude_(desiredLongitude),
            desiredRadius_(desiredRadius){ }

    ~ObjectiveFunction(){ }

    double function( ){

        double finalTime = missionSegment_->getFinalSimulationEpoch();

        std::map< double, Eigen::VectorXd > solutionMap = missionSegment_->dynamicsSimulator_
                ->getEquationsOfMotionNumericalSolution();
        std::map< double, Eigen::VectorXd >::reverse_iterator rit = solutionMap.rbegin();

        double time1 = rit->first;
        Eigen::Vector6d state1 = rit->second;

        ++rit;

        double time2 = rit->first;
        Eigen::Vector6d state2 = rit->second;

        Eigen::Vector6d finalState =  state2 + (state1 - state2)/(time1 - time2)*( finalTime - time2);

        Eigen::Vector6d finalSphericalState = orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                    finalState );

        double finalFlightPath = finalSphericalState[orbital_element_conversions::flightPathIndex];

        double finalRadius = finalSphericalState[orbital_element_conversions::radiusIndex];

        double finalLongitude = orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                    ephemerides::transformStateToTargetFrame( finalState, finalTime,
                            missionSegment_->dynamicsSimulator_->bodyMap_["Earth"]
                                    ->getRotationalEphemeris()))[orbital_element_conversions::longitudeIndex];

        //std::cout << "FinalFlightPath: " << finalFlightPath*180/mathematical_constants::PI <<"°\n";
        //std::cout << "FinalRadius: " << finalRadius/1000 << "km\nFinalLongitude: "<<
          //           finalLongitude*180/mathematical_constants::PI <<
            //         "°\n";
        //fflush(stdout);

        return fabs(finalFlightPath) + fabs(finalRadius - desiredRadius_)/100000 + fabs(finalLongitude -
                desiredLongitude_)*10;
    }

    boost::shared_ptr< optimization::MissionSegmentSettings > missionSegment_;
    double desiredLongitude_;
    double desiredRadius_;

};


BOOST_AUTO_TEST_CASE( testMissionLinker )
{

    using namespace tudat;
    using namespace tudat::optimization;
    using namespace tudat::orbital_element_conversions;

    //pagmo::random_device::set_seed(16);

    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "naif0009.tls" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    //First segment

    //Start at November 06, 2017
    double startEpoch;
    str2et_c("2017-11-06T23:00:00", &startEpoch);

    std::vector< std::string > bodiesToCreate;

    bodiesToCreate.push_back( "Earth" );

    std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings(bodiesToCreate);

    simulation_setup::NamedBodyMap bodyMap1, bodyMap2;

    Eigen::Vector6d constantEarthState = Eigen::Vector6d::Zero();

    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
            constantEarthState, "SSB", "J2000" );

    bodySettings["Earth"]->rotationModelSettings->resetOriginalFrame("J2000");

    bodyMap1 = simulation_setup::createBodies( bodySettings );

    bodyMap1["SC"] = boost::make_shared< simulation_setup::Body >();

    for( simulation_setup::NamedBodyMap::iterator it = bodyMap1.begin(); it != bodyMap1.end(); ++it)
    {
        bodyMap2[it->first] = boost::make_shared< simulation_setup::Body >(*(it->second));
    }

    double aInitial = spice_interface::getAverageRadius( "Earth" ) + 180.0E3;

    // Make a circular orbit at 180 km of altitute and i = 28.3922°
    Eigen::Vector6d initialKeplerElements = Eigen::Vector6d::Zero();
    initialKeplerElements[ semiMajorAxisIndex ] = aInitial;
    initialKeplerElements[ inclinationIndex ] = 28.3922 * mathematical_constants::PI/180;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = mathematical_constants::PI/2;

    // Create a Earth-fixed spherical coordinate system from the circular orbit
    Eigen::Vector6d initialSphericalState = convertCartesianToSphericalOrbitalState(
            ephemerides::transformStateToTargetFrame( convertKeplerianToCartesianElements( initialKeplerElements,
                    bodyMap1["Earth"]->getGravityFieldModel()->getGravitationalParameter() ), startEpoch,
                            bodyMap1[ "Earth" ]->getRotationalEphemeris() ) );

    // Set the initial coordinates above Cape Canaveral to simulate an Eastwards launch
    // (heading angle = pi/2)
    initialSphericalState[ SphericalOrbitalStateElementIndices::latitudeIndex ] =
                28.3922*mathematical_constants::PI/180;
    initialSphericalState[ SphericalOrbitalStateElementIndices::longitudeIndex ] =
            (-80.6077) * mathematical_constants::PI/180;

    // Transform into global frame
    Eigen::Vector6d initialCartesianState =  ephemerides::transformStateToGlobalFrame(
         convertSphericalOrbitalToCartesianState( initialSphericalState ), startEpoch,
                bodyMap1[ "Earth" ]->getRotationalEphemeris( ) );

   /* Eigen::Vector6d checkSphericalState = convertCartesianToSphericalOrbitalState( ephemerides::transformStateToTargetFrame(
                initialCartesianState, startEpoch, bodyMap1["Earth"]->getRotationalEphemeris() ) );
*/
    //Semimajor axis of Geostationary orbit:

    const double aGeoStationary = pow( pow( physical_constants::SIDEREAL_DAY/( 2 * mathematical_constants::PI ),
            2.0 ) * ( bodyMap1[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter() ), 1.0/3.0 ) ;

    const double aTransfer = (aGeoStationary + aInitial)/2;

    const double eTransfer = (aGeoStationary - aInitial)/(aGeoStationary + aInitial);

    Eigen::Vector6d transferKeplerElements = Eigen::Vector6d::Zero();

    transferKeplerElements[ orbital_element_conversions::semiMajorAxisIndex ] = aTransfer;
    transferKeplerElements[ orbital_element_conversions::eccentricityIndex ] = eTransfer;

    Eigen::Vector6d sphericalFromKepler = orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                orbital_element_conversions::convertKeplerianToCartesianElements(transferKeplerElements,
                                                                                 bodyMap1[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter()));

    const double transferV = sphericalFromKepler[orbital_element_conversions::speedIndex];

    const double transferTime = mathematical_constants::PI*(pow( pow( aTransfer, 3.0)/bodyMap1[ "Earth" ]->getGravityFieldModel()->getGravitationalParameter(),
                                                           0.5) );

    const double nEarth = 2*mathematical_constants::PI/physical_constants::SIDEREAL_DAY;

    const double deltaLong = (-157.8583 + 80.6077 + 360.0)*mathematical_constants::PI/180;

    const double n = pow( bodyMap1["Earth"]->getGravityFieldModel()->getGravitationalParameter()/pow( aInitial, 3.0) , 0.5);

    const double parkingTime = ( deltaLong + nEarth*transferTime - mathematical_constants::PI )
            /(n - nEarth);

    std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > >
            accelerationSettings;
    accelerationSettings["Earth"].push_back( boost::make_shared< simulation_setup::AccelerationSettings >(
                                                 basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap["SC"] = accelerationSettings;

    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back("SC");

    std::vector< std::string > centralBodies;
    centralBodies.push_back("Earth");

    basic_astrodynamics::AccelerationMap accelerationModelMap1 = simulation_setup::createAccelerationModelsMap(
                bodyMap1, accelerationMap, bodiesToPropagate, centralBodies );
    basic_astrodynamics::AccelerationMap accelerationModelMap2 = simulation_setup::createAccelerationModelsMap(
                bodyMap2, accelerationMap, bodiesToPropagate, centralBodies );

    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagator1, propagator2;

    propagator1 = boost::make_shared< propagators::TranslationalStatePropagatorSettings< > >( centralBodies,
              accelerationModelMap1, bodiesToPropagate, initialCartesianState, 1000.0 );
    propagator2 = boost::make_shared< propagators::TranslationalStatePropagatorSettings< > >( centralBodies,
              accelerationModelMap2, bodiesToPropagate, Eigen::Vector6d::Zero(), 1000.0 );

    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integrator1, integrator2;

    const double fixedStepSize = 0.5;

    integrator1 = boost::make_shared< numerical_integrators::IntegratorSettings< > >(
                numerical_integrators::rungeKutta4, startEpoch, fixedStepSize );

    integrator2 = boost::make_shared< numerical_integrators::IntegratorSettings< > >( *integrator1 );

    boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > missionPart1, missionPart2;

    missionPart1 = boost::make_shared< propagators::SingleArcDynamicsSimulator< > >( bodyMap1, integrator1, propagator1, false );

    missionPart2 = boost::make_shared< propagators::SingleArcDynamicsSimulator< > >( bodyMap2, integrator2, propagator2, false );

    boost::shared_ptr< optimization::SingleDecisionVariableSettings > decisionVariableSettings1 =
            boost::make_shared< optimization::SingleDecisionVariableSettings >(
                optimization::simulation_time_decision_variable, parkingTime - 30.0, parkingTime  );

    boost::shared_ptr< optimization::SingleDecisionVariableSettings > singleDecisionVariableSettings1 =
            boost::make_shared< optimization::SingleSphericalOrbitalElementDecisionVariableSettings >(
                orbital_elements::speed, transferV - 1.0, transferV + 1.0, "SC" );

    boost::shared_ptr< optimization::SingleDecisionVariableSettings > singleDecisionVariableSettings2 =
            boost::make_shared< optimization::SingleDecisionVariableSettings >(
                optimization::simulation_time_decision_variable, transferTime - 15.0, transferTime + 15.0 );

    std::vector< boost::shared_ptr< optimization::SingleDecisionVariableSettings > > multipleDecisionVariableSettings;

    multipleDecisionVariableSettings.push_back( singleDecisionVariableSettings1 );
    multipleDecisionVariableSettings.push_back( singleDecisionVariableSettings2 );

    boost::shared_ptr< optimization::DecisionVariableSettings > decisionVariableSettings2 = boost::make_shared< optimization::
            DecisionVariableSettings >( multipleDecisionVariableSettings );

    std::vector< boost::shared_ptr< optimization::MissionSegmentSettings > > missionSegments;
    missionSegments.push_back( boost::make_shared< optimization::MissionSegmentSettings >( missionPart1, decisionVariableSettings1 ) );
    missionSegments.push_back( boost::make_shared< optimization::MissionSegmentSettings >( missionPart2, decisionVariableSettings2 ) );

    boost::shared_ptr< ObjectiveFunction > objectiveFunction =
            boost::make_shared< ObjectiveFunction >(
                missionSegments[1], -157.8583*mathematical_constants::PI/180, aGeoStationary  );

    boost::shared_ptr< optimization::ObjectiveFunctionSettings > objectiveFunctionSettings =
           boost::make_shared< optimization::UserDefinedObjectiveFunctionSettings >(
                boost::bind( &ObjectiveFunction::function, objectiveFunction), 0.0, 1.0e-3, 20 ) ;

    missionSegments[1]->objectiveFunctionSettings_ = objectiveFunctionSettings;

    boost::shared_ptr< optimization::OptimizationSettings > optimizationSettings =
            boost::make_shared< optimization::OptimizationSettings >(optimization::global_optimization, 32, true, 1 );

    optimization::MissionLinker missionLinker( missionSegments, optimizationSettings );

    std::map< double, Eigen::VectorXd > solutionMap = missionLinker.missionSegmentSettings_[1]->dynamicsSimulator_
            ->getEquationsOfMotionNumericalSolution();

    std::map< double, Eigen::VectorXd >::reverse_iterator solution =
            solutionMap.rbegin();

    double lastTime = missionSegments[1]->getFinalSimulationEpoch();

    double finalTime = solution->first;

    Eigen::Vector6d finalCartesianState = solution->second;

    ++solution;

    double penultimateTime = solution->first;

    Eigen::Vector6d penultimateCartesianState = solution->second;


    Eigen::Vector6d lastCartesianState = penultimateCartesianState +
            (finalCartesianState - penultimateCartesianState)/( finalTime - penultimateTime )*
            (lastTime - finalTime);

    Eigen::Vector6d lastSphericalState =
            orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                lastCartesianState);

    Eigen::Vector6d lastEarthFixedSphericalState =
            orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                ephemerides::transformStateToTargetFrame( lastCartesianState, finalTime,
                    bodyMap2["Earth"]->getRotationalEphemeris( ) ) );


    std::cout << "END. The analytical parking time is: " << parkingTime << "s\nThe numerical parking time is:" <<
       missionSegments[0]->decisionVariableValues_[0] << "s\nThe analytical transfer speed is: " <<
       transferV << "m/s\nThe numerical transfer speed is: " << missionSegments[1]->decisionVariableValues_[0] <<
       "m/s\nThe analytical transfer time is: " << transferTime << "s\nThe numerical transfer time is: " <<
       missionSegments[1]->decisionVariableValues_[1] << "s\n";
    std::cout << "\nThe desired final longitude is: -157.8583°\nThe reached longitude is: "
              << lastEarthFixedSphericalState[ orbital_element_conversions::longitudeIndex ]*180/mathematical_constants::PI <<
                 "°\nThe desired radius is: " << aGeoStationary/1000 <<
                 "km\nThe reached radius is: " <<
                 lastSphericalState[orbital_element_conversions::radiusIndex]/1000 <<
                 "km\nThe final flight path angle is: " <<
                 lastSphericalState[orbital_element_conversions::flightPathIndex]*180/mathematical_constants::PI
                    <<"°\n";


}

BOOST_AUTO_TEST_SUITE_END( )

}

}
