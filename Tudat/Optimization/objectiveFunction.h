#ifndef TUDAT_OPTIMIZATION_OBJECTIVE_FUNCTION_H
#define TUDAT_OPTIMIZATION_OBJECTIVE_FUNCTION_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalstateparameters.h"

namespace tudat{

namespace optimization{

/*!
 * \brief The ObjectiveFunctionType enum
 */
enum ObjectiveFunctionType
{
    finalCartesianComponent = 0,
    finalKeplerOrbitalElement = 1,
    finalSphericalOrbitalElement = 2,
    fromFinalDependentVariable = 3,
    fromMinOrMaxDependentVariable = 4,
    userDefinedObjectiveFunction = 5
};

/*!
 * \brief The ObjectiveFunctionSettings struct
 */
struct ObjectiveFunctionSettings{

    /*!
     * \brief ObjectiveFunctionSettings
     * \param objectiveVariable
     * \param objectiveValue
     */
    ObjectiveFunctionSettings( ObjectiveFunctionType objectiveVariable, const double objectiveValue,
                               const double tolerance = 1e-3 ) :
        objectiveVariable_( objectiveVariable ), objectiveValue_( objectiveValue ), tolerance_( tolerance ),
        objectiveValueIsSet_( true ){ }

    /*!
     * \brief ObjectiveFunctionSettings
     * \param objectiveVariable
     * \param minimizeOrMaximize
     */
    ObjectiveFunctionSettings( ObjectiveFunctionType objectiveVariable, const bool minimizeOrMaximize = true ) :
        objectiveVariable_(objectiveVariable), minimizeOrMaximize_( minimizeOrMaximize ),
        objectiveValueIsSet_( false ){ }

    virtual ~ObjectiveFunction( ){ }

    ObjectiveFunctionType objectiveVariable_;
    double objectiveValue_;
    double tolerance_;
    bool minimizeOrMaximize_;

protected:

    bool objectiveValueIsSet_;

};

/*!
 * \brief The UserDefinedObjectiveFunctionSettings struct
 */
struct UserDefinedObjectiveFunctionSettings : private ObjectiveFunctionSettings{

    UserDefinedObjectiveFunctionSettings( boost::function< double > userDefinedFunction,
                                          const double objectiveValue = 0.0 , const double tolerance = 1e-3 ) :
        ObjectiveFunctionSettings( userDefinedObjectiveFunction, objectiveValue, tolerance ),
        userDefinedFunction_( userDefinedFunction ) { }

    UserDefinedObjectiveFunctionSettings( boost::function< double > userDefinedFunction,
                                          const bool minimizeOrMaximize ) :
        ObjectiveFunctionSettings( userDefinedObjectiveFunction, minimizeOrMaximize ),
        userDefinedFunction_( userDefinedFunction ) { }

    ~UserDefinedObjectiveFunctionSettings( ){ }

    boost::function< double > userDefinedFunction_;

};


struct FinalCartesianComponentObjectiveFunctionSettings : private ObjectiveFunctionSettings{

    FinalCartesianComponentObjectiveFunctionSettings( orbital_elements::CartesianElements cartesianComponent,
            const double objectiveValue, const double tolerance = 1e-3 ) :
        ObjectiveFunctionSettings( finalCartesianComponent, objectiveValue, tolerance), cartesianComponent_( cartesianComponent )
    { }

    FinalCartesianComponentObjectiveFunctionSettings( orbital_elements::CartesianElements cartesianComponent,
            const bool minimizeOrMaximize = true ) :
        ObjectiveFunctionSettings( finalCartesianComponent, minimizeOrMaximize), cartesianComponent_( cartesianComponent )
    { }

    ~FinalCartesianComponentObjectiveFunctionSettings( ){ }

    orbital_elements::CartesianElements cartesianComponent_;

};

struct FinalKeplerElementObjectiveFunctionSettings : private ObjectiveFunctionSettings{

    FinalKeplerElementObjectiveFunctionSettings( orbital_elements::KeplerianElements keplerElement,
            std::string centralBody, const double objectiveValue, const double tolerance = 1e-3 ) :
        ObjectiveFunctionSettings( finalKeplerOrbitalElement, objectiveValue, tolerance),
        keplerElement_( keplerElement ), centralBody_(centralBody)
    { }

    FinalKeplerElementObjectiveFunctionSettings( orbital_elements::KeplerianElements keplerElement,
            std::string centralBody, const bool minimizeOrMaximize = true ) :
        ObjectiveFunctionSettings( finalKeplerOrbitalElement, minimizeOrMaximize),
        keplerElement_( keplerElement ), centralBody_(centralBody)
    { }

    ~FinalKeplerElementObjectiveFunctionSettings( ){ }

    orbital_elements::KeplerianElements keplerElement_;
    std::string centralBody_;

};


struct FinalSphericalOrbitalElementObjectiveFunctionSettings : private ObjectiveFunctionSettings{

    FinalSphericalOrbitalElementObjectiveFunctionSettings( orbital_elements::SphericalOrbitalStateElements spericalOrbitalElement,
            const double objectiveValue, const double tolerance = 1e-3 ) :
        ObjectiveFunctionSettings( finalSphericalOrbitalElement, objectiveValue, tolerance), sphericalOrbitalElement_( spericalOrbitalElement )
    { }

    FinalSphericalOrbitalElementObjectiveFunctionSettings( orbital_elements::SphericalOrbitalStateElements spericalOrbitalElement,
            const bool minimizeOrMaximize = true ) :
        ObjectiveFunctionSettings( finalSphericalOrbitalElement, minimizeOrMaximize), sphericalOrbitalElement_( spericalOrbitalElement )
    { }

    ~FinalSphericalOrbitalElementObjectiveFunctionSettings( ){ }

    orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement_;

};

struct ObjectiveFunctionFromFinalDependentVariableSettings : private ObjectiveFunctionSettings{

    ObjectiveFunctionFromFinalDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const double objectiveValue, const int indexInVectorialVariable = 0, const double tolerance = 1e-3) :
    ObjectiveFunctionSettings( fromFinalDependentVariable, objectiveValue, tolerance ),
    dependentVariableSaveSettings_(dependentVariableSaveSettings), indexInVectorialVariable_(indexInVectorialVariable){ }

    ObjectiveFunctionFromFinalDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const bool minimizeOrMaximize, const int indexInVectorialVariable = 0) :
    ObjectiveFunctionSettings( fromFinalDependentVariable, minimizeOrMaximize ),
    dependentVariableSaveSettings_(dependentVariableSaveSettings), indexInVectorialVariable_(indexInVectorialVariable){ }

    virtual ~ObjectiveFunctionFromFinalDependentVariableSettings( ){ }

    int getPositionInDependentVariablesSaveVector( boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator );

    boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings_;
    int indexInVectorialVariable_;

};


struct ObjectiveFunctionFromMinOrMaxDependentVariableSettings : private ObjectiveFunctionFromFinalDependentVariableSettings{

    ObjectiveFunctionFromMinOrMaxDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const double objectiveValue, const int indexInVectorialVariable = 0, const bool minimumOrMaximum = true,
            const double tolerance = 1e-3) :
        ObjectiveFunctionFromFinalDependentVariableSettings( dependentVariableSaveSettings, objectiveValue, indexInVectorialVariable,
                                                             tolerance ),
        minimumOrMaximum_(minimumOrMaximum)
    {
        objectiveVariable_ = fromMinOrMaxDependentVariable;
    }

    ObjectiveFunctionFromMinOrMaxDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const bool minimizeOrMaximize = true, const bool minimumOrMaximum = true) :
        ObjectiveFunctionFromFinalDependentVariableSettings( dependentVariableSaveSettings, minimizeOrMaximize, indexInVectorialVariable ),
        minimumOrMaximum_( minimumOrMaximum )
     {
            objectiveVariable_ = fromMinOrMaxDependentVariable;
     }

    ~ObjectiveFunctionFromMinOrMaxDependentVariableSettings( ){ }

    double getMinOrMaxSaveVariable(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator );


    const bool minimumOrMaximum_;

};





}

}
#endif
