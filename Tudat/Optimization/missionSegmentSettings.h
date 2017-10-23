/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OPTIMIZATION_MISSION_SEGMENT_SETTINGS_H
#define TUDAT_OPTIMIZATION_MISSION_SEGMENT_SETTINGS_H

#include"decisionVariableSettings.h"
#include"objectiveFunction.h"


namespace tudat{

namespace optimization{


//! Class for the definition of a mission segment
struct MissionSegmentSettings{

    //! Constructor
    /*!
     * Constructor for the definition of multiple decision variable
     * \param dynamicsSimulator pre-defined single arc dynamics simulator that will be automatically
     * modified during the optimization
     * \param decisionVariableSettings object containing a vector of multiple decision variable settings
     * associated with the dynamics simulator
     * \param objectiveFunctionSettings object containing the objective function settings
     */
    MissionSegmentSettings( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                    boost::shared_ptr< DecisionVariableSettings > decisionVariableSettings = NULL,
                    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings = NULL) :
        dynamicsSimulator_( dynamicsSimulator ), decisionVariableSettings_( decisionVariableSettings ),
        objectiveFunctionSettings_( objectiveFunctionSettings ) { }

    //! Constructor
    /*!
     * Constructor for the definition of a single decision variable settings
     * \param dynamicsSimulator pre-defined single arc dynamics simulator that will be automatically
     * modified during the optimization
     * \param singleDecisionVariableSettings object containing a single decision variable settings
     * \param objectiveFunctionSettings object containing the objective function settings
     */
    MissionSegmentSettings( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                    boost::shared_ptr< SingleDecisionVariableSettings > singleDecisionVariableSettings,
                    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings = NULL ) :
        dynamicsSimulator_( dynamicsSimulator ), decisionVariableSettings_(
                boost::make_shared< DecisionVariableSettings >( singleDecisionVariableSettings ) ),
                objectiveFunctionSettings_( objectiveFunctionSettings ) { }

    //! Constructor
    /*!
     * Constructor for a mission segment without decision variables
     * \param dynamicsSimulator pre-defined single arc dynamics simulator that will be automatically
     * modified during the optimization
     * \param objectiveFunctionSettings object containing the objective function settings
     */
    MissionSegmentSettings( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings ) :
        dynamicsSimulator_( dynamicsSimulator ), decisionVariableSettings_( NULL ),
                objectiveFunctionSettings_( objectiveFunctionSettings ) { }

    ~MissionSegmentSettings( ){ }

    //! Storage variable for the dynamics simulator
    boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator_;

    //! Storage variable for the (multiple) decision variable settings
    boost::shared_ptr< DecisionVariableSettings > decisionVariableSettings_;

    //! Storage variable for the objective function settingd
    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings_;

    //! Storage variable for the best decision variables values
    //! after optimization
    Eigen::VectorXd decisionVariableValues_;

    //! Method to check the propagator settings in the dynamics simulator
    /*!
     * Method to check wether the dynamics simulator possesses a propagator settings with
     * state type propagatorType
     * \param propagatorType state type of the propagator settings to look for
     * \return true if the propagator with propagatorType state type is found, false otherwise
     */
    bool hasPropagatorSettings( propagators::IntegratedStateType propagatorType =
            propagators::transational_state );

    //! Method to retrieve the propagator settings from the dynamics simulator
    /*!
     * Method to retrieve the propagator settings from the dynamics simulator, specifying what
     * state type you need
     * \param propagatorType state type of the propagator settings that you need
     * \return shared pointer to the propagator settings required. Runtime error if the
     * propagator settings specified do not exist.
     */
    boost::shared_ptr< propagators::PropagatorSettings< > >  selectPropagatorSettings(
            propagators::IntegratedStateType propagatorType = propagators::transational_state);


    //! Method to retrieve the initial translational cartesian states of a certain body
    /*!
     * Method to retrieve the initial orbital states of a certain body specified by its name in the bodyMap
     * in the dynamics simulator. The method also returns the initial position of the state in the
     * initial state vector.
     * \param selectedBodyName Name of the body for which the initial orbital state is required, as
     * specified in the bodyMap keys.
     * \param initialSegmentSize Return value of the initial segment size. (specify any int)
     * \return Orbital state
     */
    Eigen::Vector6d getInitialOrbitalState( const std::string selectedBodyName, int* initialSegmentSize );


    //! Method to modify the initial orbital state of a body
    /*!
     * Method to modify the initial orbital state of a body in the translational state propagator settings
     * by specifing the position in the vector of initial states of the first member
     * \param modifiedState New orbital state to be substituted in the vector of orbital states
     * \param initialSegmentSize Position in the vector of the first member of the orbital state
     */
    void resetOrbitalState( Eigen::Vector6d modifiedState, int initialSegmentSize );


    //! Method to modify the value of a specified initial cartesian component
    /*!
     * Method to modify a specified initial cartesian component of a specified body from the bodyMap in
     * the dynamics simulator.
     * \param cartesianComponent Name of the cartesian component to be modified
     * \param selectedBodyName Name of the associated body in the bodyMap in the dynamics simulator
     * \param newValue New value of the cartesian component
     */
    void modifyCartesianComponent( orbital_elements::CartesianElements cartesianComponent,
                                   const std::string selectedBodyName, const double newValue );


    //! Method to modify the value of a specified initial Kepler element
    /*!
     * Method to modify the value a specified initial Kepler component of a specified body from the bodyMap in
     * the dynamics simulator.
     * \param keplerElement Name of the Kepler orbital element to be modified
     * \param selectedBodyName Name of the associated body in the bodyMap in the dynamics simulator
     * \param centralBody Name of the associated central body in the bodyMap in the dyamics simulator
     * \param newValue New value of the Kepler element
     */
    void modifyKeplerElement(orbital_elements::KeplerianElements keplerElement,
                             const std::string selectedBodyName, const std::string centralBody,
                             const double newValue );


    //! Method to modify the value of a specified initial spherical orbital element
    /*!
     * Method to modify the value a specified initial spherical orbital component of a specified body from the bodyMap in
     * the dynamics simulator.
     * \param sphericalOrbitalElement Name of the spherical orbital element to be modified
     * \param selectedBodyName Name of the associated body in the bodyMap in the dynamics simulator
     * \param newValue New value of the spherical orbital element
     */
    void modifySphericalComponent( orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
                                   const std::string selectedBodyName, const double newValue );

    //! Returns the final simulation epoch
    /*!
     * Returns the final simulation epoch, which is either defined in the time termination settings
     * or is interpolated from the actual termination conditions.
     * \return Final simulation epoch
     */
    double getFinalSimulationEpoch( );

    //! Returns the final interpolated orbital states
    /*!
     * Returns the final orbital states interpolated to the epoch finalTime. The orbital states
     * are stored in a map with the associated body names as key
     * \param finalTime Time to which the orbital states must be interpolated
     * \return Map of final orbital states interpolated to the epoch finalTime with the associated
     * body names as key
     */
    std::map< std::string, Eigen::Vector6d > getFinalOrbitalStates( const double finalTime );

    //! Returns the final body masses
    /*!
     * Returns the final body masses interpolated to the epoch finalTime. The masses are stored in
     * a map with the associated body names as key
     * \param finalTime Time to which the body masses must be interpolated
     * \return Map of final body masses interpolated to the epoch finalTime with the associated
     * body names as key
     */
    std::map< std::string, double > getFinalBodyMasses( const double finalTime );

    //! Reset the initial cartesian states of all the bodies
    /*!
     * Reset the initial cartesian states of all the bodies in the translational state propagator
     * settings in the dynamics simulator
     * \param newStates Map of new cartesian states with associated body name as key
     */
    void resetInitialOrbitalStates( std::map< std::string, Eigen::Vector6d > newStates );

    //! Reset the simulation start epoch in the integrator settings
    /*!
     * Reset the simulation start epoch in the integrator settings in the dynamics simulator
     * \param startEpoch New value for the simulation start epoch
     */
    void resetStartEpoch( const double startEpoch );

    //! Reset the initial body masses
    /*!
     * Reset the constant body masses and the initial body masses (if a mass propagator
     * is present) of all the bodies
     * \param newMasses Map of new body mass values with associated body names as key
     */
    void resetBodyMasses( std::map< std::string, double > newMasses );


    //! Reset simulation start epoch, initial orbital states and initial body masses
    /*!
     * Reset simulation start epoch in the integrator settings, initial orbital states in the
     * translational state propagator settings, body masses in the bodyMass and initial body
     * masses in the mass propagator settings if present
     * \param startEpoch New value of the simulation start epoch
     * \param newStates Map containing the new initial cartesian orbital states with associated body names as
     * key
     * \param newMasses Map containing the new initial body masses  with associated body names as
     * key
     */
    void setPreliminaryConditions( const double startEpoch, std::map< std::string, Eigen::Vector6d > newStates,
                                   std::map< std::string, double > newMasses = std::map< std::string, double >() );

    //! Reset the values of the decision variables in the dynamics simulator
    /*!
     * Reset the values of the decision variables in the dynamics simulator
     * \param newValues Vector of size n = number of decision variable settings
     */
    void setPropagationConditions( Eigen::VectorXd& newValues );

    //! Return the value of the objective function with the current configuration
    /*!
     * Return the value of the objective function with the current configuration
     * \return Current value of the objective function
     */
    double getObjectiveFunctionValue( void );

    Eigen::VectorXd getInitialStates( void );

};

}

}




#endif
