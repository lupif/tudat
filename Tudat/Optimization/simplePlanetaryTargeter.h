/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OPTIMIZATION_OBJECTIVE_FUNCTION_H
#define TUDAT_OPTIMIZATION_OBJECTIVE_FUNCTION_H

#include"Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h"
#include"Tudat/Optimization/missionLinker.h"

namespace tudat{

namespace optimization{

//! Class to obtain a 3-d analytical solution from the multi-revolution Lambert targeter.
struct AnalyticalSimplePlanetaryTargeter{

    //! Constructor
    /*!
     * Constructor.
     * \param departurePosition Initial coordinates [m] w.r.t. the center of gravitation
     * \param targetPosition Coordinates of target [m] w.r.t. the center of gravitation
     * \param travelTime Desired time of travel in seconds
     * \param centralBodyGravitationalParameter Grevitational parameter of the central massive body
     * \param numberOfRevolutions Desired number of orbits before rendez-vous
     * \param rightBranch Set to true for a positive transfer orbit, set to false for a retrograde
     * transfer orbit
     */
    AnalyticalSimplePlanetaryTargeter( Eigen::Vector3d departurePosition, Eigen::Vector3d targetPosition,
                             const double travelTime, const double centralBodyGravitationalParameter,
                             const int numberOfRevolutions = 0, const bool rightBranch = true ) : departurePosition_( departurePosition ),
        targetPosition_( targetPosition ), travelTime_( travelTime ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        numberOfRevolutions_(numberOfRevolutions), isRightBranch_( rightBranch )
    {

        //Methods to retrieve the orbital plane of the target orbit
        //1. Non-normalized transfer plane orthogonal vector H_unnorm = r_1 cross r_2
        Eigen::Vector3d H_unnorm = departurePosition_.cross( targetPosition_ );

        Eigen::Quaterniond quaternions;
        Eigen::Vector3d e;
        double theta;

        //2. Z-vector in original frame orientation = [0 0 1] transpose
        Eigen::Vector3d z = Eigen::Vector3d::Zero();
        z[2] = 1;

        //3. Limit on the vector size: if |H_unnorm| << 1 then the normalization
        //can suffer from double precision bias.
        //In order to avoid this see steps 4.2 - 7.2
        if( H_unnorm.norm() > 1e-9 )
        {
            //Case where |H_unnorm| > 10^-9 (arbitrary)

            //4.1. Normalized plane orthogonal vector h.
            Eigen::Vector3d h = H_unnorm.normalized();

            //5.1. Assume always positive rotation, i.e. z-component of h
            //must be positive
            if( h[2] < 0 )
                h = -h;

            //6.1. Calculate eigen vector: e = h cross (h_x h_y 0)-transposed
            Eigen::Vector3d h_xy = Eigen::Vector3d::Zero();
            h_xy.segment(0,2) = h.segment(0,2);
            e = h.cross( h_xy ).normalized();

            //7.1. Calculate rotation angle w.r.t. e
            theta = acos( h.dot( z ) );
        }

        else{

            //Case where |H_unnorm| < 10^-9: in this case the departure and target
            //points are most certainly opposite w.r.t. the center of gravitation

            //4.2. The Eigen vector lies on the x-y plane and it's perpendicular to
            //r_1 and r_2. Non normalized eigen vector: E_unnorm = (-r_1_y r_1_x 0)-transposed.
            Eigen::Vector3d E_unnorm = Eigen::Vector3d::Zero();
            E_unnorm[0] = -departurePosition_[1];
            E_unnorm[1] = departurePosition_[0];

            //5.2. The eigen vector e is normalized from E_unnorm
            e = E_unnorm.normalized();

            //6.2. r_1_xy = (r_1_x r_1_y 0)-transposed (in plane projection of departure point vector)
            Eigen::Vector3d departurePosition_xy = Eigen::Vector3d::Zero();
            departurePosition_xy.segment(0,2) = departurePosition_.segment(0,2);

            //7.2. Calculate rotation angle w.r.t. e, theta = -acos( r_1 dot-product r_1_xy )
            theta = -acos( departurePosition_.normalized().dot( departurePosition_xy.normalized() ) );
        }

        //8. Define quaternions for the definition of the rotation matrix
        quaternions.w() = cos( theta/2 );
        quaternions.x() = e[0] * sin( theta/2 );
        quaternions.y() = e[1] * sin( theta/2 );
        quaternions.z() = e[2] * sin( theta/2 );

        //9. Calculate rotation matrix from quaternions
        Eigen::Matrix3d rotationM = quaternions.normalized().toRotationMatrix();

        //10. Redefine r_1 and r_2 in the new frame orientation by multiplication with rotation matrix
        Eigen::Vector3d inPlaneDeparturePosition = rotationM.inverse() * departurePosition_;
        Eigen::Vector3d inPlaneTargetPosition = rotationM.inverse() * targetPosition_;

        //11. Nullify the (very small z~0) vertical component in the new frame orientation
        inPlaneDeparturePosition[2] = 0.0;
        inPlaneTargetPosition[2] = 0.0;

        //12. Use the MultiRevolutionLambertTargeterIzzo to calculate the in-plane velocities
        mission_segments::MultiRevolutionLambertTargeterIzzo lambertTargeter(
                    inPlaneDeparturePosition, inPlaneTargetPosition,
                    travelTime_, centralBodyGravitationalParameter_ );

        const int maxNumberOfRevolutions = lambertTargeter.getMaximumNumberOfRevolutions();

        if( numberOfRevolutions_ > maxNumberOfRevolutions )
        {
            throw std::runtime_error( "The specified number of revolutions is larger than the maximum" );
        }

        lambertTargeter.computeForRevolutionsAndBranch( numberOfRevolutions_, isRightBranch_ );

        //13. Translate velocities to original frame orientation by inverse multiplication with
        //rotation matrix.
        departureVelocity_ = rotationM * lambertTargeter.getInertialVelocityAtDeparture();
        arrivalVelocity_ = rotationM * lambertTargeter.getInertialVelocityAtArrival();

    }

    //! Cartesian departure position in original frame
    Eigen::Vector3d departurePosition_;

    //! Cartesian target position in original frame
    Eigen::Vector3d targetPosition_;

    //! Transfer time in seconds
    double travelTime_;

    //! Gravitational parameter of center of gravitation
    double centralBodyGravitationalParameter_;

    //! Number of complete orbits
    int numberOfRevolutions_;

    //! Set to true for a positive transfer orbit, set to false for a retrograde
    //! transfer orbit
    bool isRightBranch_;

    //! Departure velocity cartesian components
    Eigen::Vector3d departureVelocity_;

    //! Arrival velocity cartesian components
    Eigen::Vector3d arrivalVelocity_;

};


//! Class to numerically optimize the unperturbed transfer orbit between two planets
struct UnperturbedNumericalPlanetaryTargeter
{

    //! Constructor
    /*!
     * Constructor. This constructor will propagate and optimize the transfer between two planets
     * centered at SSB and with SPICE ephemerides of the planets and Sun, using only the central
     * gravity of the Sun as acceleration
     * \param departurePlanet Name of the departure planet in SPICE
     * \param targetPlanet Name of the target planet in SPICE
     * \param departureEpoch Departure epoch in seconds from J2000
     * \param arrivalEpoch Arrival epoch in seconds from J2000
     * \param integrationStep Integration step in seconds
     * \param referenceFrame Name of the inertial reference frame orientation
     * \param numberOfRevolutions Number of orbits before rendez-vous
     * \param rightBranch Set to true for a positive transfer orbit, set to false for a retrograde
     * transfer orbit
     */
    UnperturbedNumericalPlanetaryTargeter( std::string departurePlanet, std::string targetPlanet,
            const double departureEpoch, const double arrivalEpoch, double integrationStep,
            std::string referenceFrame = "ECLIPJ2000", const int numberOfRevolutions = 0, const bool rightBranch = true )
    {

        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( departurePlanet );
        bodiesToCreate.push_back( targetPlanet );


        std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettings =
                simulation_setup::getDefaultBodySettings( bodiesToCreate, departureEpoch - 300, arrivalEpoch + 300 );

        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( referenceFrame );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( referenceFrame );
        }

        simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "SC" ] = boost::make_shared< simulation_setup::Body >( );

        // Finalize body creation.
        simulation_setup::setGlobalFrameBodyEphemerides( bodyMap, "SSB", referenceFrame );

        // Define propagator settings variables.
        simulation_setup::SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bodiesToPropagate.push_back( "SC" );
        centralBodies.push_back( "SSB" );

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfSC;
        accelerationsOfSC[ "Sun" ].push_back( boost::make_shared< simulation_setup::AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationMap[ "SC" ] = accelerationsOfSC;

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        Eigen::Vector6d departureState = bodyMap[ departurePlanet ]->
                getStateInBaseFrameFromEphemeris( departureEpoch );

        Eigen::Vector6d arrivalState = bodyMap[ targetPlanet ]->
                getStateInBaseFrameFromEphemeris( arrivalEpoch );

        Eigen::Vector3d sunCenteredDeparture = departureState.segment(0, 3) -
                bodyMap[ "Sun" ]->getStateInBaseFrameFromEphemeris( departureEpoch ).segment(0, 3);

        Eigen::Vector3d sunCenteredArrival = arrivalState.segment(0, 3) -
                bodyMap[ "Sun" ]->getStateInBaseFrameFromEphemeris( departureEpoch ).segment(0, 3);

        //Eigen::Vector3d sunVelocityInFrame = bodyMap[ "Sun" ]->getStateInBaseFrameFromEphemeris( departureEpoch ).segment(3,3);

        AnalyticalSimplePlanetaryTargeter analyticalSolution( sunCenteredDeparture, sunCenteredArrival, arrivalEpoch - departureEpoch,
                bodyMap[ "Sun" ]->getGravityFieldModel()->getGravitationalParameter(), numberOfRevolutions, rightBranch );

        departureState.segment(3,3) = analyticalSolution.departureVelocity_;// + sunVelocityInFrame;

        // Create numerical integrator.
        boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
                boost::make_shared< numerical_integrators::IntegratorSettings< > >
                ( numerical_integrators::rungeKutta4, departureEpoch, integrationStep );

        boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > distanceFromTarget =
                boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::relative_distance_dependent_variable, "SC", targetPlanet );

        std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back( distanceFromTarget );

        boost::shared_ptr< propagators::DependentVariableSaveSettings > depVarSaveSettings=
                boost::make_shared< propagators::DependentVariableSaveSettings >( dependentVariables, false);

        boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, departureState, arrivalEpoch,
                  propagators::cowell, depVarSaveSettings );

        Eigen::VectorXd lowerBoundaries;

        lowerBoundaries.resize(3);
        lowerBoundaries[0] = departureState[3] - 100;
        lowerBoundaries[1] = departureState[4] - 100;
        lowerBoundaries[2] = departureState[5] - 100;

        Eigen::VectorXd upperBoundaries;

        upperBoundaries.resize(3);
        upperBoundaries[0] = departureState[3] + 100;
        upperBoundaries[1] = departureState[4] + 100;
        upperBoundaries[2] = departureState[5] + 100;

        boost::shared_ptr< optimization::SingleDecisionVariableSettings > decisionVarSettings =
                boost::make_shared< optimization::SingleDecisionVariableSettings >(
                    optimization::initial_cartesian_velocity_decision_variable,
                    lowerBoundaries, upperBoundaries, "SC" );

        boost::shared_ptr< optimization::ObjectiveFunctionFromFinalDependentVariableSettings > objectiveFuntionSettings
                = boost::make_shared< optimization::ObjectiveFunctionFromFinalDependentVariableSettings >(
                        distanceFromTarget, 0, 0, spice_interface::getAverageRadius( targetPlanet ), 10 );

        // Create simulation object and propagate dynamics.
        boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator =
                boost::make_shared< propagators::SingleArcDynamicsSimulator< > >(
                    bodyMap, integratorSettings, propagatorSettings, false );

        boost::shared_ptr< MissionSegmentSettings > missionSegmentSettings =
                boost::make_shared< MissionSegmentSettings >( dynamicsSimulator, decisionVarSettings,
                                                              objectiveFuntionSettings );
        std::vector< boost::shared_ptr< MissionSegmentSettings > > vecOfMissionSegments;
        vecOfMissionSegments.push_back( missionSegmentSettings );

        MissionLinker missionLinker( vecOfMissionSegments,
                                     boost::make_shared< OptimizationSettings >( global_optimization, 64, true, 1 ) );

        int x = 0;
        Eigen::Vector6d initialOrbitalState = missionSegmentSettings->getInitialOrbitalState( "SC", &x );

        lowerBoundaries[0] = initialOrbitalState[3] - 5;
        lowerBoundaries[1] = initialOrbitalState[4] - 5;
        lowerBoundaries[2] = initialOrbitalState[5] - 5;

        upperBoundaries[0] = initialOrbitalState[3] + 5;
        upperBoundaries[1] = initialOrbitalState[4] + 5;
        upperBoundaries[2] = initialOrbitalState[5] + 5;

        decisionVarSettings->boundaries_->setLowerBoundary( lowerBoundaries );
        decisionVarSettings->boundaries_->setUpperBoundary( lowerBoundaries );

        objectiveFuntionSettings->maxNumberOfEvolutions_ = 10;

        missionLinker = MissionLinker(vecOfMissionSegments,
                                      boost::make_shared< OptimizationSettings >( global_optimization, 64, true, 1 ) );


        initialOrbitalState = missionSegmentSettings->getInitialOrbitalState( "SC", &x );

        lowerBoundaries[0] = initialOrbitalState[3] - 0.01;
        lowerBoundaries[1] = initialOrbitalState[4] - 0.01;
        lowerBoundaries[2] = initialOrbitalState[5] - 0.01;

        upperBoundaries[0] = initialOrbitalState[3] + 0.01;
        upperBoundaries[1] = initialOrbitalState[4] + 0.01;
        upperBoundaries[2] = initialOrbitalState[5] + 0.01;

        decisionVarSettings->boundaries_->setLowerBoundary( lowerBoundaries );
        decisionVarSettings->boundaries_->setUpperBoundary( lowerBoundaries );

        objectiveFuntionSettings->maxNumberOfEvolutions_ = 10;

        missionLinker = MissionLinker(vecOfMissionSegments,
                                      boost::make_shared< OptimizationSettings >( global_optimization, 64, true, 1 ) );


        initialOrbitalState = missionSegmentSettings->getInitialOrbitalState( "SC", &x );

        Eigen::Vector6d finalOrbitalState = missionSegmentSettings->getFinalOrbitalStates( arrivalEpoch )["SC"];

        inertialDepartureVelocity_ = initialOrbitalState.segment(3,3);
        inertialArrivalVelocity_ = finalOrbitalState.segment(3,3);

        inertialDepartureDeltaV_ = inertialDepartureVelocity_ - departureState.segment(3,3);
        inertialArrivalDeltaV_ = inertialArrivalVelocity_ - arrivalState.segment(3,3);

        radialDepartureVelocity_ = inertialDepartureVelocity_.norm();
        radialArrivalVelocity_ = inertialArrivalVelocity_.norm();

        departureDeltaV_ = inertialDepartureDeltaV_.norm();
        arrivalDeltaV_ = inertialArrivalDeltaV_.norm();

        spacecraftStateMap_ = dynamicsSimulator->getEquationsOfMotionNumericalSolution();

        for( double time = departureEpoch; time < arrivalEpoch + integrationStep;
             time += integrationStep ){

            departurePlanetStateMap_[ time ] = bodyMap[ departurePlanet ]->
                    getStateInBaseFrameFromEphemeris( time );

            targetPlanetStateMap_[ time ] = bodyMap[ targetPlanet ]->
                    getStateInBaseFrameFromEphemeris( time );

        }

    }


    //! Inertial departure velocity
    Eigen::Vector3d inertialDepartureVelocity_;
    //! Inertial arrival velocity
    Eigen::Vector3d inertialArrivalVelocity_;
    //! Velocity difference at departure with departure planet
    Eigen::Vector3d inertialDepartureDeltaV_;
    //! Velocity difference at arrival with target planet
    Eigen::Vector3d inertialArrivalDeltaV_;
    //! Departure speed in SSB centered reference frame
    double radialDepartureVelocity_;
    //! Arrivl speed in SSB centered reference frame
    double radialArrivalVelocity_;
    //! Departure delta-v (w.r.t. departure planet orbit)
    double departureDeltaV_;
    //! Arrival delta-v (w.r.t. target planet orbit)
    double arrivalDeltaV_;

    //! State map of the spacecraft during transfer
    std::map< double, Eigen::VectorXd > spacecraftStateMap_;

    //! State map of departure planet during transfer
    std::map< double, Eigen::VectorXd > departurePlanetStateMap_;

    //! State map of target planet during transfer
    std::map< double, Eigen::VectorXd > targetPlanetStateMap_;

};




}

}





#endif
