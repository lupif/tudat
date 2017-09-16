#ifndef TUDAT_OPTIMIZATION_DECISION_VAR_SETTINGS

#define TUDAT_OPTIMIZATION_DECISION_VAR_SETTINGS



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
    ~Boundaries( ){ };

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
    double getLowerBoundary( const unsigned int index = 0 );

    //! Method to access a position in the multivariate upper boundary
    /*!
     * Method to access a position in the multivariate upper boundary and
     * retrieve the value. It can be used to access directly 1 sized (scalar)
     * upper boundaries
     * \param index position in the multivariate upper boundary
     * \return double value of the position index
     */
    double getUpperBoundary( const unsigned int index = 0 );


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

    simulationTime = 0,//Tunes the total simulation time (final epoch - initial epoch)
    fromTerminationSettings = 1,
    initialStateCartesianComponents = 2, //Reset completely the initial state
    initialVelocitytateCartesianComponents = 3, //Modifies the initial velocity
    initialPositionStateCartesianComponents = 4, //Modifies the initial position
    singleCartesianComponentsElement = 5, //Change any from the state
    singleKeplerOrbitalElement = 6, //Change any from kepler orbital elements
    singleSphericalOrbitalElement = 7,//Change any from spherical orbital elements

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
     * \param decisionVariable type of decision variable
     * \param boundaries boundaries
     */
    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    boost::shared_ptr< Boundaries > boundaries );

    //! Constructor
    /*!
     * Constructor accepting scalar lower and upper boundaries
     * \param decisionVariable decision variable type
     * \param lowerBoundary lower boundary for scalar decision variable
     * \param upperBoundary upper boundary for scalar decision variable
     */
    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    double lowerBoundary, double upperBoundary );

    //! Constructor
    /*!
     * Constructor accepting vectorial lower and upper boundaries
     * \param decisionVariable decision variable type
     * \param lowerBoundary lower boundary for vectorial decision variable
     * \param upperBoundary upper boundary for vectorial decision variable
     */
    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    Eigen::VectorXd& lowerBoundary, Eigen::VectorXd& upperBoundary);

    //! Deconstructor
    virtual ~SingleDecisionVariableSettings( ){ }


    DecisionVariables decisionVariable_;
    boost::shared_ptr< Boundaries > boundaries_;

};


//! Class to use one of the termination settings as decision variable (including the final time)
struct SingleDecisionVariableFromTerminationSettings : public SingleDecisionVariableSettings
{
    //! Constructor
    /*!
     * \Constructor.
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
        SingleDecisionVariableSettings( fromTerminationSettings, lowerBoundary, upperBoundary ),
        positionInVectorOfTerminationSettings_( positionInVectorOfTerminationSettings )
    {
        checkTerminationConditions(dynamicsSimulator);
    }

    ~SingleDecisionVariableFromTerminationSettings(){ }

    //! index in the vector of termination settings in
    //! case of multiple termination condition settings.
    int positionInVectorOfTerminationSettings_;

    //! Pointer used to directly address the position of the value to modify in the termination settings
    double* memoryPositionOfVariable_;


private:

    //! Scouts the termination conditions and sets the memoryPositionOfVariable_ variable
    void checkTerminationConditions(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator);

};


//! Class to define a cartesian component as decision variable
struct SingleCartesianComponentDecisionVariableSettings : public SingleDecisionVariableSettings{

    SingleCartesianComponentDecisionVariableSettings( orbital_elements::CartesianElements cartesianComponent,
                                                 double lowerBoundary, double upperBoundary ) :
        SingleDecisionVariableSettings( singleCartesianComponentsElement, lowerBoundary, upperBoundary ),
        cartesianComponent_( cartesianComponent ) { }

    ~SingleCartesianComponentDecisionVariableSettings( ){ }

    orbital_elements::CartesianElements cartesianComponent_;

};

//! Class to define a kepler element as decision variable
struct SingleKeplerElementDecisionVariableSettings : public SingleDecisionVariableSettings{

    SingleKeplerElementDecisionVariableSettings( orbital_elements::KeplerianElements keplerElement,
                                                 std::string centralBody, double lowerBoundary,
                                                 double upperBoundary ) :
        SingleDecisionVariableSettings( singleKeplerOrbitalElement, lowerBoundary, upperBoundary),
        centralBody_( centralBody ), keplerElement_(keplerElement) { }

    ~SingleKeplerElementDecisionVariableSettings(){ }

    std::string centralBody_;
    orbital_elements::KeplerianElements keplerElement_;

};

//! Class to define a spherical orbital component as decision variable
struct SingleSphericalOrbitalElementDecisionVariableSettings : public SingleDecisionVariableSettings{

    SingleSphericalOrbitalElementDecisionVariableSettings( orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
                                                 double lowerBoundary, double upperBoundary ) :
        SingleDecisionVariableSettings( singleSphericalOrbitalElement, lowerBoundary, upperBoundary),
        sphericalOrbitalElement_(sphericalOrbitalElement) { }

    ~SingleSphericalOrbitalElementDecisionVariableSettings( ){ }

    orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement_;

};


struct DecisionVariableSettings{

    DecisionVariableSettings( std::vector< boost::shared_ptr< SingleDecisionVariableSettings > >
            multiDecisionVariableSettings ) : decisionVariableSettings_( multiDecisionVariableSettings )
    { }

    std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > decisionVariableSettings_;

};


}

}

#endif
