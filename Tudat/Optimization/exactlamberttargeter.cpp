#include "exactlamberttargeter.h"



tudat::optimization::ExactLambertTargeter::ExactLambertTargeter( tudat::simulation_setup::NamedBodyMap bodyMap,
        std::string departurePlanetName, std::string targetPlanetName, std::string spacecraftName, const double departureDay,                                                         
        const double timeOfFlight, tudat::simulation_setup::SelectedAccelerationMap accelerationMap,
        tudat::propagators::TranslationalPropagatorType userDefinedPropagatorType,
        boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings > userDefinedIntegratorSettings,
        std::vector< std::string >& bodiesToPropagate, std::vector< std::string >& centralBodies,
        Eigen::VectorXd& initialStates) : bodyMap_(bodyMap), spacecraftName_(spacecraftName), timeOfFlight_(timeOfFlight),
        departureDay_(departureDay)
{

        //Scout the bodiesToPropagate to find out whether the two planets, the spacecraft and the Sun are
        //there
        Eigen::Vector6d targetPlanetCartesianComponents;
        Eigen::Vector6d sunCartesianComponents;
        bool foundSpacecraftInBodiesToPropagate = false;
        bool propagateSun = false;
        unsigned int i = 0;
        for( std::vector< std::string >::iterator it = bodiesToPropagate.begin(); it++; it != bodiesToPropagate.end() )
        {
            if( std::strcmp( *it, departurePlanetName_ ) == 0 )
            {
                departurePlanetCartesianComponents_ = initialStates.segment( (int)i * 6, 6 );
                propagateDeparturePlanet_ = true;
                departurePlanetNamePos_ = i;
            }
            if( std::strcmp( *it, targetPlanetName_ ) == 0 )
            {
                targetPlanetCartesianComponents = initialStates.segment( (int)i * 6, 6 );
                propagateTargetPlanet_ = true;
                targetPlanetNamePos_ = i;
            }
            if( std::strcmp( *it, spacecraftName_ ) == 0 )
            {
                phase1InitialState_ = initialStates.segment( (int)i * 6, 6 );
                foundSpacecraftInBodiesToPropagate = true;
                spacecraftNamePos_ = i;
            }
            if( std::strcmp( *it, "Sun" ) == 0 )
            {
                sunCartesianComponents = initialStates.segment( (int)i * 6, 6 );
                propagateSun = true;
            }

            i++;
        }

        // If Sun is not to propagate retrieve ephemeris from bodyMap
        if( !propagateSun )
        {
            sunCartesianComponents = bodyMap_[ "Sun" ]
                    ->getStateInBaseFrameFromEphemeris( departureDay );
        }

        // If departure planet is not to propagate retrieve ephemeris from bodyMap
        if( !propagateDeparturePlanet )
        {
        departurePlanetCartesianComponents_ = bodyMap_[ departurePlanetName ]
                ->getStateInBaseFrameFromEphemeris( departureDay ) - sunCartesianComponents;
        }

        // If target planet is not to propagate retrieve ephemeris from bodyMap
        if( !propagateTargetPlanet )
        {
        targetPlanetCartesianComponents = bodyMap_[ targetPlanetName ]
                ->getStateInBaseFrameFromEphemeris( departureDay ) - sunCartesianComponents;
        }

        // If spacecraft is not found in bodies to propagate set initial ephemeris at target planet
        if( !foundSpacecraftBodiesToPropagate )
        {
            bodiesToPropagate_.push_back( spacecraftName_ );
            centralBodies_.push_back( "Sun" );
            initialStates.resize( i * 6 + 6 );
            initialStates.segment( i * 6, 6 ) = departurePlanetCartesianComponents_;
        }

        /*-------------Preliminary calculations------------------*/


        const double departurePlanetGravitationalParameter = bodyMap_[ departurePlanetName ]
                ->getGravityFieldModel( )->getGravitationalParameter( );
        const double targetPlanetGravitationalParameter = bodyMap_[ targetPlanetName ]
                ->getGravityFieldModel( )->getGravitationalParameter( );
        const double sunGravitationalParameter = bodyMap_[ "Sun" ]
                ->getGravityFieldModel( )->getGravitationalParameter( );


        const Eigen::Vector6d departurePlanetKeplerianComponents =
                tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                    departurePlanetCartesianComponents, sunGravitationalParameter );

        const Eigen::Vector6d targetPlanetKeplerianComponents =
                tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                    targetPlanetCartesianComponents, sunGravitationalParameter );

        const Eigen::Vector6d spacecraftKeplerianComponents =
                tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                    phase1InitialState_, sunGravitationalParameter );

        const double departurePlanetSemimajorAxis = departurePlanetKeplerianComponents.at( tudat::orbital_element_conversions::semiMajorAxisIndex );
        const double departurePlanetEccentricity = departurePlanetKeplerianComponents.at( tudat::orbital_element_conversions::eccentricityIndex );
        const double targetPlanetSemimajorAxis = targetPlanetKeplerianComponents.at( tudat::orbital_element_conversions::semiMajorAxisIndex );
        const double targetPlanetEccentricity = targetPlanetKeplerianComponents.at( tudat::orbital_element_conversions::eccentricityIndex );
        const double spacecraftSemimajorAxis = spacecraftKeplerianComponents.at( tudat::orbital_element_conversions::semiMajorAxisIndex );
        const double spacecraftEccentricity = spacecraftKeplerianComponents.at( tudat::orbital_element_conversions::eccentricityIndex );


        departurePlanetHillSphereRadius_ =
                ( 1 - departurePlanetEccentricity ) * departurePlanetSemimajorAxis
                * std::pow( departurePlanetGravitationalParameter/( 3 * sunGravitationalParameter), 1/3 );
        targetPlanetHillSphereRadius_ =
                ( 1 - targetPlanetEccentricity ) * targetPlanetSemimajorAxis
                * std::pow( targetPlanetGravitationalParameter/( 3 * sunGravitationalParameter), 1/3 );
        maximumSimulationTime_ = std::pow( std::pow( ( ( 1 + spacecraftEccentricity ) * spacecraftSemimajorAxis +
                ( 1 + targetPlanetEccentricity ) * targetPlanetSemimajorAxis ), 3) / sunGravitationalParameter, 0.5 ) *
                tudat::mathematical_constants::PI * 1.1;


        /*-----------------Preliminary Acceleration maps---------------*/

        //Switch off all the acceleration components on the spacecraft except the Sun during passage
        //in the spheres of influence of the planets
        phase1and3AccelerationMap_ = phase2AccelerationMap_;
        tudat::simulation_setup::SelectedAccelerationMap::iterator it = phase1and3AccelerationMap_.find( spacecraftName_ );
        if( it != phase1and3AccelerationMap_.end() )
        {
            phase1and3AccelerationMap_.erase( it );
        }
        std::map< std::vector< boost::shared_ptr< tudat::simulation_setup::AccelerationSettings > > > accelerationOfSpacecraft;
        accelerationOfSpacecraft["Sun"].push_back( boost::make_shared< tudat::simulation_setup::AccelerationSettings > (
                    tudat::basic_astrodynamics::central_gravity ) );
        phase1and3AccelerationMap_[spacecraftName_] = accelerationOfSpacecraft;

        tudat::basic_astrodynamics::AccelerationMap accelerationModelMapPhase1and3 = tudat::simulation_setup::createAccelerationModelsMap(
                    bodyMap_, phase1and3AccelerationMap_, bodiesToPropagate_, centralBodies_ );
        tudat::basic_astrodynamics::AccelerationMap accelerationModelMapPhase2_ = tudat::simulation_setup::createAccelerationModelsMap(
                    bodyMap_, phase2AccelerationMap_, bodiesToPropagate_, centralBodies_ );

        /*---------------------Setup termination settings for the 3 phases-------------------*/

        //Add the distances from the planets to the dependent variables to save
        dependentVariablesList_.push_back(
                    boost::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >( tudat::propagators::relative_distance_dependent_variable, spacecraftName_, departurePlanetName_ ) );
        dependentVariablesList_.push_back(
                    boost::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >( tudat::propagators::relative_distance_dependent_variable, spacecraftName_, targetPlanetName_ ) );
        dependentVariablesToSave_ =
                        boost::make_shared< tudat::propagators::DependentVariableSaveSettings >( dependentVariablesList );


        //Distance from planets
        boost::shared_ptr< tudat::propagators::SingleDependentVariableSaveSettings > terminationDistancePhase1 =
                boost::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    tudat::propagators::relative_distance_dependent_variable, spacecraftName, departurePlanetName );
        boost::shared_ptr< tudat::propagators::SingleDependentVariableSaveSettings > terminationDistancePhase2and3 =
                boost::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    tudat::propagators::relative_distance_dependent_variable, spacecraftName, targetPlanetName );

        terminationSettingsPhase1_ =
                boost::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                        terminationDistancePhase1, departurePlanetHillSphereRadius_, false );
        terminationSettingsPhase2_.push_back(
                boost::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                        terminationDistancePhase2and3, targetPlanetHillSphereRadius_, true ) );
        terminationSettingsPhase3_ =
                boost::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                        terminationDistancePhase2and3, targetPlanetHillSphereRadius_, false );

        // propagation time
        terminationSettingsPhase2_.push_back(
                boost::make_shared< tudat::propagators::PropagationTimeTerminationSettings >(
                    maximumSimulationTime_ + simulationStartEpoch_ ) );

        // Hybrid termination settings for phase 2
        hybridTerminationSettingsPhase2_ =
                boost::make_shared< tudat::propagators::PropagationHybridTerminationSettings > ( terminationSettingsPhase2, true );


        // Integrator Settings for Phase 1
        double phase1FixedStepSize = departurePlanetHillSphereRadius_/1e9*3600;
        integratorSettingsPhase1_ =
                boost::make_shared< tudat::numerical_integrators::IntegratorSettings< > >
        ( rungeKutta4, departureDay_, phase1FixedStepSize );

}



double tudat::optimization::ExactLambertTargeter::propagatePlanetaryTransfer( Eigen::Vector3d& vInfiniteVectorial )
{


    //Initialize initial cartesian components adding additional velocity
    Eigen::Vector6d phase1InitialState= phase1InitialState_;
    phase1InitialState.segment(3,4) += vInfiniteVectorial;

    // Create propagation settings for phase 1
    // Set up a cowell propagator for phase 1
    boost::shared_ptr< tudat::propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsPhase1 =
            boost::make_shared< tudat::propagators::TranslationalStatePropagatorSettings< double > >
            ( centralBodies_, phase1and3AccelerationMap_, bodiesToPropagate_, phase1InitialState,
              terminationSettingsPhase1_, cowell, dependentVariablesToSave_ );

    tudat::propagators::SingleArcDynamicsSimulator phase1Simulator( bodyMap_,
           integratorSettingsPhase1_, propagatorSettingsPhase1, true, false, false );


    Eigen::Vector6d phase2InitialState = phase1Simulator.getEquationsOfMotionNumericalSolution().rbegin()->second;
    double phase2StartEpoch = phase1Simulator.getEquationsOfMotionNumericalSolution().rbegin()->first;

    boost::shared_ptr< tudat::propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsPhase2 =
            boost::make_shared< tudat::propagators::TranslationalStatePropagatorSettings< double > >
            ( centralBodies_, phase2AccelerationMap_, bodiesToPropagate_, systemInitialStatePhase2,
              hybridTerminationSettingsPhase2_, userDefinedPropagator_, dependentVariablesToSave_, );

    // Continue from here

    userDefinedIntegratorSettings_->initialTime_ = phase2StartEpoch;


    tudat::propagators::SingleArcDynamicsSimulator phase2Simulator( bodyMap_,
            userDefinedIntegratorSettings_, phase2PropagatorSettings, true, false, false );

    double phase3StartEpoch = phase2Simulator.getEquationsOfMotionNumericalSolution().rbegin()->first;

    if( phase3StartEpoch >= maximumSimulationTime_ )
    {
        return phase2Simulator.getDependentVariableHistory().rbegin()->second.at(0);
    }

    Eigen::Vector6d phase3InitialState = phase2Simulator.getEquationsOfMotionNumericalSolution().rbegin()->second;


    boost::shared_ptr< tudat::propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsPhase3 =
            boost::make_shared< tudat::propagators::TranslationalStatePropagatorSettings< double > >
            ( centralBodies_, phase1and3AccelerationMap_, bodiesToPropagate_, phase3InitialState,
            terminationSettingsPhase1_, cowell, dependentVariablesToSave_ );

    double phase3FixedStepSize = targetPlanetHillSphereRadius_/1e9*3600;

    boost::shared_ptr< tudat::numerical_integrators::IntegratorSettings< > > integratorSettingsPhase3 =
            boost::make_shared< tudat::numerical_integrators::IntegratorSettings< > >( rungeKutta4,
                    phase3StartEpoch, phase3FixedStepSize );

    phase2Simulator.getDependentVariableHistory().rbegin()->second;

    tudat::propagators::SingleArcDynamicsSimulator phase3Simulator( bodyMap_,
                integratorSettingsPhase3, propagatorSettingsPhase3, true, false, false );

    std::map< double, std::vector< double > > dependentVariableHistory =
            phase3Simulator.getDependentVariableHistory();

    double bestDistanceFromTarget = dependentVariableHistory.rbegin()->second.at(0);

    double distanceFromTarget;

    for( std::map< double, std::vector< double > >::reverse_iterator it = dependentVariableHistory.rbegin();
         it++; it!=dependentVariableHistory.rend() )
    {
        distanceFromTarget = it->second.at(0);
        if( distanceFromTarget < bestDistanceFromTarget )
            bestDistanceFromTarget = distanceFromTarget;
    }

    return bestDistanceFromTarget;
    
}



