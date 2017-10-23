/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OPTIMIZATION_DECISION_VAR_SETTINGS_H

#define TUDAT_OPTIMIZATION_DECISION_VAR_SETTINGS_H



#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalstateparameters.h"

namespace tudat{

namespace optimization{

//! Class used to store the upper and lower (multivariate) boundaries of the decision variables
struct Boundaries{

    //! Constructor
    /*!
     * Constructor.
     * \param lowerBoundary multivariate lower boundary
     * \param upperBoundary multivariate upper boundary
    */
    Boundaries( const Eigen::VectorXd lowerBoundary,
                const Eigen::VectorXd upperBoundary );

    //! Deconstructor
    ~Boundaries( ){ }

    //! Method to retrieve the lower boundary
    /*!
     * Get the stored lower boundary
     * \return lower boundary
     */
    Eigen::VectorXd getLowerBoundary( );

    //! Method to retrieve the upper boundary
    /*!
     * Get the stored upper boundary
     * \return upper boundary
     */
    Eigen::VectorXd getUpperBoundary( );

    //! Method to access a position in the multivariate lower boundary
    /*!
     * Method to access a position in the multivariate lower boundary and
     * retrieve the value. It can be used to access directly 1 sized (scalar)
     * lower boundaries
     * \param index position in the multivariate lower boundary
     * \return double value of the position index
     */
    double getLowerBoundary( const unsigned int index );

    //! Method to access a position in the multivariate upper boundary
    /*!
     * Method to access a position in the multivariate upper boundary and
     * retrieve the value. It can be used to access directly 1 sized (scalar)
     * upper boundaries
     * \param index position in the multivariate upper boundary
     * \return double value of the position index
     */
    double getUpperBoundary( const unsigned int index );


    //! Method to set or modify the lower boundary
    /*!
     * Method to set or modify the lower boundary.
     * \param lowerBoundary new (multivariate) lower boundary
     */
    void setLowerBoundary( const Eigen::VectorXd lowerBoundary );


    //! Method to set or modify the upper boundary
    /*!
     * Method to set or modify the upper boundary.
     * \param upperBoundary new (multivariate) upper boundary
     */
    void setUpperBoundary( const Eigen::VectorXd upperBoundary );


private:

    //! Store variable of lower boundary
    Eigen::VectorXd lowerBoundary_;

    //! Store variable of upper boundary
    Eigen::VectorXd upperBoundary_;

};



//! Enum of decision variables types
//! \brief The enum of decision variables types.
enum DecisionVariables
{

    simulation_time_decision_variable = 0, //Tunes the total simulation time (final epoch - initial epoch)
    from_termination_settings_decision_variable = 1, //Use this for e.g. final simulation time and other termination settings
    initial_cartesian_state_decision_variable = 2, //Reset completely the initial state
    initial_cartesian_velocity_decision_variable = 3, //Modifies the initial velocity
    initial_cartesian_position_decision_variable = 4, //Modifies the initial position
    single_cartesian_component_decision_variable = 5, //Change any from the state
    single_kepler_element_decision_variable = 6, //Change any from kepler orbital elements
    single_spherical_orbital_element_decision_variable = 7,//Change any from spherical orbital elements

};


//! Function to retrieve the size of the decision variable by name
/*!
 * Function to retrieve the size of the decision variable by name
 * \param decisionVariable enumerator of the decision variable
 * \return int size of the decision variable (1 for scalar >1 if vector,
 * -1 if undetermined)
 */
int getDecisionVariableSize( DecisionVariables decisionVariable );



//! Class to define the settings for a decision variable
struct SingleDecisionVariableSettings{

    //! Constructor
    /*!
     * Constructor accepting boundaries as class Boundaries.
     * \param decisionVariable Name of decision variable
     * \param boundaries Boundaries object
     * \param associatedBody Name of body to which the decision variable
     * belongs.
     */
    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    boost::shared_ptr< Boundaries > boundaries,
                                    std::string associatedBody = "" );

    //! Constructor
    /*!
     * Constructor accepting scalar lower and upper boundaries
     * \param decisionVariable decision variable type
     * \param lowerBoundary lower boundary for scalar decision variable
     * \param upperBoundary upper boundary for scalar decision variable
     * \param associatedBody Name of body to which the decision variable
     * belongs.
     */
    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    double lowerBoundary, double upperBoundary,
                                    std::string associatedBody = "" );

    //! Constructor
    /*!
     * Constructor accepting vectorial lower and upper boundaries
     * \param decisionVariable decision variable type
     * \param lowerBoundary lower boundary for vectorial decision variable
     * \param upperBoundary upper boundary for vectorial decision variable
     * \param associatedBody Name of body to which the decision variable
     * belongs.
     */
    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    Eigen::VectorXd& lowerBoundary, Eigen::VectorXd& upperBoundary,
                                    std::string associatedBody = "" );

    //! Deconstructor
    virtual ~SingleDecisionVariableSettings( ){ }

    //! Decision variable name
    DecisionVariables decisionVariable_;

    //! Box-boundaries object
    boost::shared_ptr< Boundaries > boundaries_;

    //! Name of the body for which the decision variable holds
    std::string associatedBody_;
};

//! Class to use one of the termination settings as decision variable (including the final time)
struct SingleDecisionVariableFromTerminationSettings : public SingleDecisionVariableSettings
{
    //! Constructor
    /*!
     * Constructor.
     * \param dynamicsSimulator dynamics simulator in which the termination settings belong.
     * \param lowerBoundary lower boundary of the termination setting.
     * \param upperBoundary upper boundary of the termination setting.
     * \param positionInVectorOfTerminationSettings index in the vector of termination settings in
     * case of multiple termination condition settings.
     */
    SingleDecisionVariableFromTerminationSettings(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
            const double lowerBoundary, const double upperBoundary,
            int positionInVectorOfTerminationSettings = 0 ) :
        SingleDecisionVariableSettings( from_termination_settings_decision_variable, lowerBoundary, upperBoundary ),
        positionInVectorOfTerminationSettings_( positionInVectorOfTerminationSettings )
    {
        checkTerminationConditions(dynamicsSimulator);
    }

    ~SingleDecisionVariableFromTerminationSettings(){ }

    //! Index in the vector of termination settings in
    //! case of multiple termination condition settings.
    int positionInVectorOfTerminationSettings_;

    //! Pointer used to directly address the position of the value to modify in the termination settings
    double* memoryPositionOfVariable_;


private:

    //! Scouts the termination conditions in a SingleArcDynamicsSimulator
    //! and sets the memoryPositionOfVariable_ variable
    /*!
     * Scouts the termination conditions and sets the memoryPositionOfVariable_ variable
     * \param dynamicsSimulator SingleArcDyanmicsSimulator obejct in which the termination
     * conditions to be checked are stored.
     */
    void checkTerminationConditions(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator);

};


//! Class to define a cartesian component as decision variable
struct SingleCartesianComponentDecisionVariableSettings : public SingleDecisionVariableSettings{

    //! Constructor
    /*!
     * Constructor.
     * \param cartesianComponent Cartesian component chosen among the six defining the state
     * \param lowerBoundary Lower boundary of the decision variable
     * \param upperBoundary Upper boundary of the decision variable
     * \param associatedBody String containing the name of the simulated body as stored in the
     * dynamics simulator's bodyMap.
     */
    SingleCartesianComponentDecisionVariableSettings( orbital_elements::CartesianElements cartesianComponent,
                                                 double lowerBoundary, double upperBoundary, std::string associatedBody ) :
        SingleDecisionVariableSettings( single_cartesian_component_decision_variable, lowerBoundary, upperBoundary,
                                        associatedBody ),
        cartesianComponent_( cartesianComponent ) { }

    //! Deconstructor
    ~SingleCartesianComponentDecisionVariableSettings( ){ }

    //! Cartesian component chosen among the six defining the state
    orbital_elements::CartesianElements cartesianComponent_;

};

//! Class to define a kepler element as decision variable
struct SingleKeplerElementDecisionVariableSettings : public SingleDecisionVariableSettings{

    //! Constructor
    /*!
     * Constructor.
     * \param keplerElement Kepler element chosen among the 6 defining orbit and anomaly
     * \param centralBody name of the central massive body as stored in the dynamics simulator's bodyMap
     * around which the (osculating) orbit is defined
     * \param lowerBoundary lower boundary of decision variable
     * \param upperBoundary upper boundary of decision variable
     * \param associatedBody name of the simulated body as stored in the
     * dynamics simulator's bodyMap.
     */
    SingleKeplerElementDecisionVariableSettings( orbital_elements::KeplerianElements keplerElement,
                                                 std::string centralBody, double lowerBoundary,
                                                 double upperBoundary, std::string associatedBody ) :
        SingleDecisionVariableSettings( single_kepler_element_decision_variable, lowerBoundary, upperBoundary, associatedBody),
        centralBody_( centralBody ), keplerElement_(keplerElement) { }

    ~SingleKeplerElementDecisionVariableSettings(){ }

    //! centralBody name of the central massive body as stored in the dynamics simulator's bodyMap
    //! around which the (osculating) orbit is defined
    std::string centralBody_;

    //!  Kepler element chosen among the 6 defining orbit and anomaly
    orbital_elements::KeplerianElements keplerElement_;

};

//! Class to define a spherical orbital component as decision variable
struct SingleSphericalOrbitalElementDecisionVariableSettings : public SingleDecisionVariableSettings{

    //! Constructor
    /*!
     * Constructor.
     * \param sphericalOrbitalElement Spherical orbital element chosen among the 6 defining the spherical
     * orbital state
     * \param lowerBoundary Lower boundary of decision variable
     * \param upperBoundary Upper boundary of decision variable
     * \param associatedBody Name of the simulated body as stored in the
     * dynamics simulator's bodyMap.
     */
    SingleSphericalOrbitalElementDecisionVariableSettings( orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
                                                 double lowerBoundary, double upperBoundary, std::string associatedBody ) :
        SingleDecisionVariableSettings( single_spherical_orbital_element_decision_variable, lowerBoundary, upperBoundary, associatedBody),
        sphericalOrbitalElement_(sphericalOrbitalElement) { }

    //! Deconstructor
    ~SingleSphericalOrbitalElementDecisionVariableSettings( ){ }

    //! Name of spherical orbital element chosen among the 6 defining the spherical
    //! orbital state.
    orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement_;

};


//! Class containing multiple SingleDecisionVariableSettings
struct DecisionVariableSettings{

    //! Constructor
    /*!
     * Constructor accepting a vector of SingleDecisionVariableSettings objects
     * \param multiDecisionVariableSettings vector of single decision variable settings
     */
    DecisionVariableSettings( std::vector< boost::shared_ptr< SingleDecisionVariableSettings > >
            multiDecisionVariableSettings ) : decisionVariableSettings_( multiDecisionVariableSettings )
    { }

    //! Constructor
    /*!
     * Constructor accepting a SingleDecisionVariableSettings object
     * \param singleDecisionVariableSettings Single decision variable settings object
     */
    DecisionVariableSettings( boost::shared_ptr< SingleDecisionVariableSettings >
            singleDecisionVariableSettings )
    {
        decisionVariableSettings_.push_back( singleDecisionVariableSettings );
    }

    //! Vector of SingleDecisionVariableSettings.
    std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > decisionVariableSettings_;

};


}

}

#endif
