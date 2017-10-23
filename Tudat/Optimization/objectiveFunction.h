/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OPTIMIZATION_OBJECTIVE_FUNCTION
#define TUDAT_OPTIMIZATION_OBJECTIVE_FUNCTION

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalstateparameters.h"
#include "boost/function.hpp"
#include "boost/bind.hpp"

namespace tudat{

namespace optimization{

/*!
 * \brief The ObjectiveFunctionType enum
 */
enum ObjectiveFunctionType
{
    final_cartesian_component_objective_function = 0,
    final_kepler_orbital_element_objective_function = 1,
    final_spherical_orbital_element_objective_function = 2,
    final_dependent_variable_objective_function = 3,
    min_max_dependent_variable_objective_function = 4,
    user_defined_objective_function = 5
};

//!  Base class defining the objective function of the optimization process
struct ObjectiveFunctionSettings{

    //! Constructor
    /*!
     * Constructor for the definition of an objective function with objective convergence value.
     * \param objectiveVariable Name of the variable to be used as objective function
     * \param objectiveValue User defined convergence value
     * \param tolerance Error from the objective value for which the optimization stops evolving
     * \param maxNumberOfEvolutions Maximum number of iterations in the optimization
     */
    ObjectiveFunctionSettings( ObjectiveFunctionType objectiveVariable,
                               const double objectiveValue, const double tolerance = 1e-3,
                               const unsigned int maxNumberOfEvolutions = 100 ) :
        objectiveVariable_( objectiveVariable ), objectiveValue_( objectiveValue ), tolerance_( tolerance ),
        maxNumberOfEvolutions_( maxNumberOfEvolutions ), objectiveValueIsSet_( true ){ }

    //! Constructor
    /*!
     * Constructor for the definition of an objective function with options to minimize or maximize continuosly
     * the objective function without any specified convergence value
     * \param objectiveVariable Name of the variable to be used as objective function
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize
     * \param maxNumberOfEvolutions Maximum number of iterations in the optimization
     */
    ObjectiveFunctionSettings( ObjectiveFunctionType objectiveVariable, const bool minimizeOrMaximize = true,
                               const unsigned int maxNumberOfEvolutions = 100 ) :
        objectiveVariable_(objectiveVariable), minimizeOrMaximize_( minimizeOrMaximize ),
        maxNumberOfEvolutions_( maxNumberOfEvolutions ), objectiveValueIsSet_( false ){ }

    virtual ~ObjectiveFunctionSettings( ){ }

    //! Name of the variable to be used as objective function
    ObjectiveFunctionType objectiveVariable_;

    //! User defined convergence value
    double objectiveValue_;

    //! Error from the objective value for which the optimization stops evolving
    double tolerance_;

    //! Variable set to true to minimize, set to false to maximize the objective function
    bool minimizeOrMaximize_;

    //! Maximum number of iterations in the optimization
    unsigned int maxNumberOfEvolutions_;

    //! This variable is set to True if the first contructor is used
    bool objectiveValueIsSet_;


};


//! Class to define as objective function the final value of one of the 6 Cartesian components
struct FinalCartesianComponentObjectiveFunctionSettings : public ObjectiveFunctionSettings{

    //! Constructor
    /*!
     * Constructor with option to set an objective value
     * \param cartesianComponent Name of the cartesian component ( [x,y,z]CartesianPosition or [x,y,z]CartesianVelocity )
     * \param associatedBody Name of the body for which the cartesian component is retrieved, as stored in the bodyMap
     * of the dynamics simulator
     * \param objectiveValue Desired convergence value for the objective function
     * \param tolerance Error for which the convergence is achieved
     * \param maxNumberOfEvolutions Maximum number of optimization iterations
     */
    FinalCartesianComponentObjectiveFunctionSettings( orbital_elements::CartesianElements cartesianComponent,
            const std::string associatedBody, const double objectiveValue, const double tolerance = 1e-3,
                                                      const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( final_cartesian_component_objective_function, objectiveValue, tolerance,
                maxNumberOfEvolutions ), cartesianComponent_( cartesianComponent ),
        associatedBody_(associatedBody)
    { }


    //! Constructor
    /*!
     * Constructor used to define whether to minimize or maximize the function.
     * \param cartesianComponent Name of the cartesian component ( [x,y,z]CartesianPosition or [x,y,z]CartesianVelocity )
     * \param associatedBody Name of the body for which the cartesian component is retrieved, as stored in the bodyMap of
     * the dynamics simulator.
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize
     * \param maxNumberOfEvolutions Maximum number of optimization iterations
     */
    FinalCartesianComponentObjectiveFunctionSettings( orbital_elements::CartesianElements cartesianComponent,
            const std::string associatedBody, const bool minimizeOrMaximize = true, const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( final_cartesian_component_objective_function, minimizeOrMaximize,
                maxNumberOfEvolutions ), cartesianComponent_( cartesianComponent ),
        associatedBody_(associatedBody)
    { }

    ~FinalCartesianComponentObjectiveFunctionSettings( ){ }

    //! Name of the cartesian component
    orbital_elements::CartesianElements cartesianComponent_;

    //! Name of the body for which the cartesian component is retrieved
    std::string associatedBody_;

};

//! Class to define the final value of one of the 6 Kepler elements of a body as objective function
struct FinalKeplerElementObjectiveFunctionSettings : public ObjectiveFunctionSettings{


    //! Constructor
    /*!
     * Constructor with option to set an objective value.
     * \param keplerElement Name of the Kepler element.
     * \param associatedBody Name of the body for which the Kepler element is retrieved, as stored in the bodyMap
     * of the dynamics simulator
     * \param centralBody Name of the central body around which the kepler elements are computed, as stored in the bodyMap of
     * the dynamics simulator.
     * \param objectiveValue Desired convergence value for the objective function.
     * \param tolerance Error for which the convergence is achieved.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    FinalKeplerElementObjectiveFunctionSettings( orbital_elements::KeplerianElements keplerElement,
            const std::string associatedBody, std::string centralBody, const double objectiveValue,
            const double tolerance = 1e-3, const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( final_kepler_orbital_element_objective_function, objectiveValue,
                tolerance, maxNumberOfEvolutions ),
        keplerElement_( keplerElement ),
        associatedBody_(associatedBody), centralBody_(centralBody)
    { }


    //! Constructor
    /*!
     * Constructor used to define whether to minimize or maximize the function.
     * \param keplerElement Name of the kepler element
     * \param associatedBody Name of the body for which the Kepler element is retrieved, as stored in the bodyMap of
     * the dynamics simulator.
     * \param centralBody Name of the central body around which the kepler elements are computed, as stored in the bodyMap of
     * the dynamics simulator.
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    FinalKeplerElementObjectiveFunctionSettings( orbital_elements::KeplerianElements keplerElement,
            const std::string associatedBody, std::string centralBody, const bool minimizeOrMaximize = true,
            const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( final_kepler_orbital_element_objective_function,
                minimizeOrMaximize, maxNumberOfEvolutions ),
        keplerElement_( keplerElement ),
        associatedBody_(associatedBody), centralBody_(centralBody)
    { }

    ~FinalKeplerElementObjectiveFunctionSettings( ){ }

    //! Name of the kepler element.
    orbital_elements::KeplerianElements keplerElement_;

    //! Name of the body for which the Kepler element is retrieved.
    std::string associatedBody_;

    //! Name of the central body around which the kepler elements are computed.
    std::string centralBody_;

};


//! Class to define the final value of one of the 6 spherical orbital elements of a body as objective function
struct FinalSphericalOrbitalElementObjectiveFunctionSettings : public ObjectiveFunctionSettings{

    //! Constructor
    /*!
     * Constructor with option to set an objective value.
     * \param spericalOrbitalElement Name of the spherical orbital element
     * \param associatedBody Name of the body for which the spherical orbital element is retrieved, as stored in the bodyMap of
     * the dynamics simulator.
     * \param objectiveValue Desired convergence value for the objective function.
     * \param tolerance Error for which the convergence is achieved.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    FinalSphericalOrbitalElementObjectiveFunctionSettings( orbital_elements::SphericalOrbitalStateElements spericalOrbitalElement,
            const std::string associatedBody, const double objectiveValue, const double tolerance = 1e-3,
            const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( final_spherical_orbital_element_objective_function,
                objectiveValue, tolerance, maxNumberOfEvolutions ), sphericalOrbitalElement_( spericalOrbitalElement ),
        associatedBody_(associatedBody)
    { }

    //! Constructor
    /*!
     * Constructor used to define whether to minimize or maximize the function.
     * \param spericalOrbitalElement Name of the spherical orbital element.
     * \param associatedBody Name of the body for which the spherical orbital element is retrieved, as stored in the bodyMap of
     * the dynamics simulator.
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    FinalSphericalOrbitalElementObjectiveFunctionSettings( orbital_elements::SphericalOrbitalStateElements spericalOrbitalElement,
            const std::string associatedBody, const bool minimizeOrMaximize = true,
            const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( final_spherical_orbital_element_objective_function,
                minimizeOrMaximize, maxNumberOfEvolutions ), sphericalOrbitalElement_( spericalOrbitalElement ),
        associatedBody_(associatedBody)
    { }

    ~FinalSphericalOrbitalElementObjectiveFunctionSettings( ){ }

    //! Name of the spherical orbital element.
    orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement_;

    //! Name of the body for which the spherical orbital element is retrieved.
    std::string associatedBody_;
};


//! Class to define the final value of one of the dependent variables as the objective function
struct ObjectiveFunctionFromFinalDependentVariableSettings : public ObjectiveFunctionSettings{

    //! Constructor
    /*!
     *  Constructor with option to set an objective value.
     * \param dependentVariableSaveSettings SingleDependentVariableSettings object as stored in the dependentVariableSaveSettings in the
     * propagator whose final value is to be used as objective function.
     * \param objectiveValue Desired convergence value for the objective function.
     * \param indexInVectorialVariable Index of the desired value in a vectorial or multivariate dependent variable.
     * \param tolerance Error for which the convergence is achieved.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    ObjectiveFunctionFromFinalDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const double objectiveValue, const int indexInVectorialVariable = 0, const double tolerance = 1e-3,
            const unsigned int maxNumberOfEvolutions = 100 ) :
    ObjectiveFunctionSettings( final_dependent_variable_objective_function, objectiveValue, tolerance,
            maxNumberOfEvolutions ),
    dependentVariableSaveSettings_(dependentVariableSaveSettings), indexInVectorialVariable_(indexInVectorialVariable){ }

    //! Costructor
    /*!
     * Constructor used to define whether to minimize or maximize the function.
     * \param dependentVariableSaveSettings SingleDependentVariableSettings object as stored in the dependentVariableSaveSettings in the
     * propagator whose final value is to be used as objective function.
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize.
     * \param indexInVectorialVariable Index of the desired value in a vectorial or multivariate dependent variable.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    ObjectiveFunctionFromFinalDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const bool minimizeOrMaximize, const int indexInVectorialVariable = 0,
            const unsigned int maxNumberOfEvolutions = 100 ) :
    ObjectiveFunctionSettings( final_dependent_variable_objective_function, minimizeOrMaximize,
            maxNumberOfEvolutions ),
    dependentVariableSaveSettings_(dependentVariableSaveSettings), indexInVectorialVariable_(indexInVectorialVariable){ }

    virtual ~ObjectiveFunctionFromFinalDependentVariableSettings( ){ }

    //! Method to retrieve the position in the vector of dependent variables to save of the desired value
    /*!
     * Retrieve the position of the desired value in the vector of dependent variables to save
     * \param dynamicsSimulator SingleArcDynamicsSimulator object in which the dependent variables to save have been defined
     * \return Index in the vector of dependent variables to save if a match with the dependentVariableSaveSettings_
     * member has been found, -1 otherwise.
     */
    int getPositionInDependentVariablesSaveVector( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator );

    //! SingleDependentVariableSettings object as stored in the dependentVariableSaveSettings in the
    //! propagator whose final value is to be used as objective function.
    boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings_;

    //! Index of the desired value in a vectorial or multivariate dependent variable.
    int indexInVectorialVariable_;

};


//! Class to define the minimum or maximum value of one of the dependent variables in the history of dependent variables to save as the objective function
struct ObjectiveFunctionFromMinOrMaxDependentVariableSettings : public ObjectiveFunctionFromFinalDependentVariableSettings{

    //! Constructor
    /*!
     * Constructor with option to set an objective value.
     * \param dependentVariableSaveSettings SingleDependentVariableSettings object as stored in the dependentVariableSaveSettings in the
     * propagator whose final value is to be used as objective function.
     * \param objectiveValue Desired convergence value for the objective function.
     * \param indexInVectorialVariable Index of the desired value in a vectorial or multivariate dependent variable.
     * \param minimumOrMaximum Set to true to retrieve the minimum value, set to false to retrieve the maximum value
     * \param tolerance Error for which the convergence is achieved.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    ObjectiveFunctionFromMinOrMaxDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const double objectiveValue, const int indexInVectorialVariable = 0, const bool minimumOrMaximum = true,
            const double tolerance = 1e-3, const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionFromFinalDependentVariableSettings( dependentVariableSaveSettings, objectiveValue, indexInVectorialVariable,
                tolerance, maxNumberOfEvolutions ),
        minimumOrMaximum_(minimumOrMaximum)
    {
        objectiveVariable_ = min_max_dependent_variable_objective_function;
    }

    //! Constructor
    /*!
     *
     * Constructor used to define whether to minimize or maximize the function.
     * \param dependentVariableSaveSettings SingleDependentVariableSettings object as stored in the dependentVariableSaveSettings in the
     * propagator whose final value is to be used as objective function.
     * \param indexInVectorialVariable Index of the desired value in a vectorial or multivariate dependent variable.
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize.
     * \param minimumOrMaximum Set to true to retrieve the minimum value, set to false to retrieve the maximum value
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    ObjectiveFunctionFromMinOrMaxDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const int indexInVectorialVariable = 0, const bool minimizeOrMaximize = true, const bool minimumOrMaximum = true,
            const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionFromFinalDependentVariableSettings( dependentVariableSaveSettings,
                minimizeOrMaximize, indexInVectorialVariable, maxNumberOfEvolutions ),
        minimumOrMaximum_( minimumOrMaximum )
     {
            objectiveVariable_ = min_max_dependent_variable_objective_function;
     }

    ~ObjectiveFunctionFromMinOrMaxDependentVariableSettings( ){ }

    //! Retrieve the minimum or maximum dependent variable defined by the user in the constructor
    //! from the history map (run the propagation before using this function)
    /*!
     * \param dynamicsSimulator SingleArcDynamicsSimulator object with the dependent variable
     * history map.
     * \return value of the minimum or maximum (as defined by the user in the constructor) of the
     * variable defined by the user in the constructor in the dependent variable history map
     */
    double getMinOrMaxSaveVariable(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator );


    //! This member is set to true to retrieve the minimum value, set to false to retrieve the maximum value
    const bool minimumOrMaximum_;

};

//! Class for defining a custom function as objective function
struct UserDefinedObjectiveFunctionSettings : public ObjectiveFunctionSettings{

    //! Constructor
    /*!
     * Constructor with option to define an objective value
     * \param userDefinedFunction Double function with no parameters used to calculate the obejctive function.
     * \param objectiveValue Desired convergence value for the objective function.
     * \param tolerance Error for which the convergence is achieved.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    UserDefinedObjectiveFunctionSettings( boost::function< double () > userDefinedFunction,
                                          const double objectiveValue = 0.0 , const double tolerance = 1e-3,
                                          const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( user_defined_objective_function, objectiveValue, tolerance,
                                   maxNumberOfEvolutions ),
        userDefinedFunction_( userDefinedFunction ) { }


    //! Constructor
    /*!
     * Constructor used to define whether to minimize or maximize the function.
     * \param userDefinedFunction Double function with no parameters used to calculate the obejctive function.
     * \param minimizeOrMaximize Set to true to minimize, set to false to maximize.
     * \param maxNumberOfEvolutions Maximum number of optimization iterations.
     */
    UserDefinedObjectiveFunctionSettings( boost::function< double () > userDefinedFunction,
                                          const bool minimizeOrMaximize,
                                          const unsigned int maxNumberOfEvolutions = 100 ) :
        ObjectiveFunctionSettings( user_defined_objective_function, minimizeOrMaximize,
                                   maxNumberOfEvolutions ),
        userDefinedFunction_( userDefinedFunction ) { }

    ~UserDefinedObjectiveFunctionSettings( ){ }

    //!Double function with no parameters used to calculate the obejctive function.
    boost::function< double () > userDefinedFunction_;

};


//! Compare two SingleDependentVariableSettings objects
/*!
 * \param singleDependentVariableSaveSettings1 SingleDependentVariableSettings first object
 * \param singleDependentVariableSaveSettings2 SingleDependentVariableSettings second object
 * \return true if the two objects have been defined with the same constructor, false otherwise
 */
bool compareSingleDependentVariableSettings(
        boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings1,
        boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings2 );

}

}

#endif
