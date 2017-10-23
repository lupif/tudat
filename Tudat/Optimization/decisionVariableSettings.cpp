/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include"decisionVariableSettings.h"

namespace tudat{

namespace optimization{

Boundaries::Boundaries( const Eigen::VectorXd lowerBoundary,
                        const Eigen::VectorXd upperBoundary ) :
    lowerBoundary_(lowerBoundary), upperBoundary_(upperBoundary) { }

Eigen::VectorXd Boundaries::getLowerBoundary(){

    return lowerBoundary_;
}

Eigen::VectorXd Boundaries::getUpperBoundary(){

    return upperBoundary_;
}

double Boundaries::getLowerBoundary( const unsigned int index ){

    return lowerBoundary_[ index ];
}

double Boundaries::getUpperBoundary( const unsigned int index ){

    return upperBoundary_[ index ];
}

void Boundaries::setLowerBoundary( const Eigen::VectorXd lowerBoundary ){

    lowerBoundary_ = lowerBoundary;

}

void Boundaries::setUpperBoundary( const Eigen::VectorXd upperBoundary ){

    upperBoundary_ = upperBoundary;
}


int getDecisionVariableSize( DecisionVariables decisionVariable )
{

    int size = -1;
    switch( decisionVariable ){
    case simulation_time_decision_variable:
        size = 1;
        break;
    case initial_cartesian_state_decision_variable:
        size = 6;
        break;
    case initial_cartesian_velocity_decision_variable:
        size = 3;
        break;
    case initial_cartesian_position_decision_variable:
        size = 3;
        break;
    case single_cartesian_component_decision_variable:
        size = 1;
        break;
    case single_kepler_element_decision_variable:
        size = 1;
        break;
    case single_spherical_orbital_element_decision_variable:
        size = 1;
        break;
    case from_termination_settings_decision_variable:
        break;
    default:
        std::string errorMessage = "Error, did not recognize decision variable size of type: " +
                boost::lexical_cast< std::string >( decisionVariable );
        throw std::runtime_error( errorMessage );
    }

    return size;
}



SingleDecisionVariableSettings::SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                boost::shared_ptr< Boundaries > boundaries,
                                std::string associatedBody ) :
    decisionVariable_(decisionVariable), boundaries_( boundaries ), associatedBody_( associatedBody ){ }


SingleDecisionVariableSettings::SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                double lowerBoundary, double upperBoundary, std::string associatedBody ) :
    decisionVariable_(decisionVariable), associatedBody_( associatedBody )
{
    if( getDecisionVariableSize( decisionVariable ) != 1 && decisionVariable != from_termination_settings_decision_variable ){
        std::string message = boost::lexical_cast< std::string >( decisionVariable ) +
                "is a vector of " + boost::lexical_cast< std::string >( getDecisionVariableSize( decisionVariable ) ) +
                ". Lower boundary and upper Boundary should be vectors.";
        throw std::runtime_error( message );
    }
    else{
        //Create boundaries from scalar input boundaries
        Eigen::VectorXd upperBoundaryVector, lowerBoundaryVector;
        upperBoundaryVector.resize(1);
        upperBoundaryVector[0] = upperBoundary;
        lowerBoundaryVector.resize(1);
        lowerBoundaryVector[0] = lowerBoundary;
        boundaries_ = boost::make_shared< Boundaries >( lowerBoundaryVector, upperBoundaryVector );
    }
}

SingleDecisionVariableSettings::SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                Eigen::VectorXd& lowerBoundary, Eigen::VectorXd& upperBoundary,
                                std::string associatedBody) :
    decisionVariable_(decisionVariable), associatedBody_( associatedBody )
{
    if( lowerBoundary.rows() != upperBoundary.rows() ){
        throw std::runtime_error( "Lower boundary and upper boundary size mismatch" );
    }
    else if( getDecisionVariableSize( decisionVariable ) != lowerBoundary.rows()
             && decisionVariable != from_termination_settings_decision_variable ){
        std::string message = "Mismatch: decision variable " + boost::lexical_cast< std::string >( decisionVariable ) +
                " has " + boost::lexical_cast< std::string >( getDecisionVariableSize( decisionVariable ) ) +
                " components, while defined boundaries have " + boost::lexical_cast< std::string >( lowerBoundary.rows() ) + ".";
        throw std::runtime_error ( message );
    }
    else
        boundaries_ = boost::make_shared< Boundaries >( lowerBoundary, upperBoundary );
}



void SingleDecisionVariableFromTerminationSettings::checkTerminationConditions(
        boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator)
{

    propagators::PropagationTerminationTypes terminationType;

    boost::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings =
        dynamicsSimulator->getPropagatorSettings()->getTerminationSettings();

    if( terminationSettings->terminationType_ == propagators::time_stopping_condition ){

        boost::shared_ptr< propagators::PropagationTimeTerminationSettings > timeTerminationSettings =
                boost::dynamic_pointer_cast< propagators::PropagationTimeTerminationSettings >(
                dynamicsSimulator->getPropagatorSettings()->getTerminationSettings() );

        memoryPositionOfVariable_ = &(timeTerminationSettings->terminationTime_);
    }
    else if( terminationSettings->terminationType_ == propagators::dependent_variable_stopping_condition ){
        boost::shared_ptr< propagators::PropagationDependentVariableTerminationSettings >
                depVarTerminationSettings =
                boost::dynamic_pointer_cast< propagators::PropagationDependentVariableTerminationSettings >(
                    dynamicsSimulator->getPropagatorSettings()->getTerminationSettings() );
        memoryPositionOfVariable_ = &(depVarTerminationSettings->limitValue_);

    }
    else if( terminationSettings->terminationType_ == propagators::hybrid_stopping_condition ){
        boost::shared_ptr< propagators::PropagationHybridTerminationSettings > hybridTerminationSettings =
                 boost::dynamic_pointer_cast< propagators::PropagationHybridTerminationSettings >(
                    dynamicsSimulator->getPropagatorSettings()->getTerminationSettings());
        terminationType = hybridTerminationSettings->
                terminationSettings_[ positionInVectorOfTerminationSettings_ ]->terminationType_;
        if( terminationType == propagators::dependent_variable_stopping_condition )
        {
            boost::shared_ptr< propagators::PropagationDependentVariableTerminationSettings >
                    depVarTerminationSettings = boost::dynamic_pointer_cast<
                    propagators::PropagationDependentVariableTerminationSettings >( hybridTerminationSettings->
                    terminationSettings_[positionInVectorOfTerminationSettings_] );
            memoryPositionOfVariable_ = &(depVarTerminationSettings->limitValue_);
        }
        else if( terminationType == propagators::time_stopping_condition )
        {
            boost::shared_ptr< propagators::PropagationTimeTerminationSettings > timeTerminationSettings =
                    boost::dynamic_pointer_cast< propagators::PropagationTimeTerminationSettings > (
                        hybridTerminationSettings->terminationSettings_[positionInVectorOfTerminationSettings_] );
            memoryPositionOfVariable_ = &( timeTerminationSettings->terminationTime_ );
        }
    }
}


} //

}
