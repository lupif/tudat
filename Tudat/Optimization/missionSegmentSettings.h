#ifndef TUDAT_OPTIMIZATION_MISSION_SEGMENT_SETTINGS_H
#define TUDAT_OPTIMIZATION_MISSION_SEGMENT_SETTINGS_H

#include"decisionVariableSettings.h"
#include"objectiveFunction.h"


namespace tudat{

namespace optimization{


struct MissionSegmentSettings{


    MissionSegmentSettings( boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator,
                    boost::shared_ptr< DecisionVariableSettings > decisionVariableSettings = NULL,
                    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings = NULL) :
        dynamicsSimulator_( dynamicsSimulator ), decisionVariableSettings_( decisionVariableSettings ),
        objectiveFunctionSettings_( objectiveFunctionSettings ) { }

    MissionSegmentSettings( boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator,
                    boost::shared_ptr< SingleDecisionVariableSettings > singleDecisionVariableSettings = NULL,
                    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings = NULL ) :
        dynamicsSimulator_( dynamicsSimulator ), decisionVariableSettings_(
                boost::make_shared< DecisionVariableSettings >( {singleDecisionVariableSettings} ) ),
                objectiveFunctionSettings_( objectiveFunctionSettings ) { }

    ~MissionSegmentSettings( ){ }

    boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator_;
    boost::shared_ptr< DecisionVariableSettings > decisionVariableSettings_;
    boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings_;

    Eigen::Vector6d getInitialOrbitalState( const std::string selectedBodyName, int* initialSegmentSize );

    void resetOrbitalState( Eigen::Vector6d& modifiedState, int initialSegmentSize );

    void modifyCartesianComponent( orbital_elements::CartesianElements cartesianComponent,
                                   const std::string selectedBodyName, const double newValue );


    void modifyKeplerElement(orbital_elements::KeplerianElements keplerElement,
                             const std::string selectedBodyName, const std::string centralBody,
                             const double newValue );


    void modifySphericalComponent( orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
                                   const std::string selectedBodyName, const double newValue );

    double getFinalSimulationEpoch( );


};








}
}




#endif
