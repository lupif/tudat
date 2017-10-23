/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include"objectiveFunction.h"

namespace tudat{

namespace optimization{


bool compareSingleDependentVariableSettings(
        boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings1,
        boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings2 )
{

    bool isEqual = false;

    if( ( singleDependentVariableSaveSettings1->variableType_ ==
          singleDependentVariableSaveSettings2->variableType_ ) &&
                ( singleDependentVariableSaveSettings1->associatedBody_ ==
                singleDependentVariableSaveSettings2->associatedBody_ ) &&
                        ( singleDependentVariableSaveSettings1->secondaryBody_ ==
                         singleDependentVariableSaveSettings2->secondaryBody_ ) )
    {
        if( ( singleDependentVariableSaveSettings1->variableType_ ==
                propagators::single_acceleration_norm_dependent_variable ) ||
                        ( singleDependentVariableSaveSettings1->variableType_ ==
                                propagators::single_acceleration_dependent_variable ) )
        {
            boost::shared_ptr< propagators::SingleAccelerationDependentVariableSaveSettings >
                    accelerationVarSaveSettings1 =
                            boost::dynamic_pointer_cast<
                                    propagators::SingleAccelerationDependentVariableSaveSettings >(
                                            singleDependentVariableSaveSettings1 );
            boost::shared_ptr< propagators::SingleAccelerationDependentVariableSaveSettings >
                    accelerationVarSaveSettings2 =
                            boost::dynamic_pointer_cast<
                                    propagators::SingleAccelerationDependentVariableSaveSettings >(
                                            singleDependentVariableSaveSettings2 );
            if( accelerationVarSaveSettings1->accelerationModeType_ ==
                    accelerationVarSaveSettings2->accelerationModeType_)
                isEqual = true;
        }
        else if( singleDependentVariableSaveSettings1->variableType_ ==
                propagators::intermediate_aerodynamic_rotation_matrix_variable )
        {
            boost::shared_ptr< propagators::IntermediateAerodynamicRotationVariableSaveSettings >
                    intermAerodRotVarSaveSettings1 =
                            boost::dynamic_pointer_cast<
                                    propagators::IntermediateAerodynamicRotationVariableSaveSettings >(
                                            singleDependentVariableSaveSettings1 );
            boost::shared_ptr< propagators::IntermediateAerodynamicRotationVariableSaveSettings >
                    intermAerodRotVarSaveSettings2 =
                            boost::dynamic_pointer_cast<
                                    propagators::IntermediateAerodynamicRotationVariableSaveSettings >(
                                            singleDependentVariableSaveSettings2);
            if( ( intermAerodRotVarSaveSettings1->baseFrame_ == intermAerodRotVarSaveSettings2->baseFrame_)
                    && (intermAerodRotVarSaveSettings1->targetFrame_ ==
                            intermAerodRotVarSaveSettings2->targetFrame_ ) )
                isEqual = true;
        }
        else if( singleDependentVariableSaveSettings1->variableType_ ==
                propagators::relative_body_aerodynamic_orientation_angle_variable )
        {
            boost::shared_ptr< propagators::BodyAerodynamicAngleVariableSaveSettings >
                    bodyAerodAngleVarSaveSettings1 =
                            boost::dynamic_pointer_cast<
                                    propagators::BodyAerodynamicAngleVariableSaveSettings >(
                                            singleDependentVariableSaveSettings1 );
            boost::shared_ptr< propagators::BodyAerodynamicAngleVariableSaveSettings >
                    bodyAerodAngleVarSaveSettings2 =
                            boost::dynamic_pointer_cast<
                                    propagators::BodyAerodynamicAngleVariableSaveSettings >(
                                            singleDependentVariableSaveSettings2);
            if( bodyAerodAngleVarSaveSettings1->angle_ == bodyAerodAngleVarSaveSettings2->angle_ )
                isEqual = true;
        }
        else
            isEqual = true;
    }

    return isEqual;

}


int ObjectiveFunctionFromFinalDependentVariableSettings::getPositionInDependentVariablesSaveVector(
        boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator ){

    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
                dynamicsSimulator->getPropagatorSettings() );

    propagators::PropagationDependentVariables desiredDependentVariableType =
            dependentVariableSaveSettings_->variableType_;

    if( indexInVectorialVariable_ >= propagators::getDependentVariableSize( desiredDependentVariableType  ) )
        throw std::runtime_error( "The desired position in the chosen variable is higher than its size" );

    int sizeOfDepVarsToSaveVector = propagatorSettings->getDependentVariablesToSave()->
            dependentVariables_.size();
    int i = 0;
    int positionInVector = 0;
    bool found = false;
    do{
        propagators::PropagationDependentVariables currentDependentVariableType =
                propagatorSettings->getDependentVariablesToSave()->
                            dependentVariables_[i]->variableType_;
        found = compareSingleDependentVariableSettings( dependentVariableSaveSettings_,
                                                        propagatorSettings->getDependentVariablesToSave()->
                                                                   dependentVariables_[i]);
        if( !found  )
            positionInVector += propagators::getDependentVariableSize( currentDependentVariableType );

        i++;
    }
    while( found == false && i < sizeOfDepVarsToSaveVector );

    if( found == false && i >= sizeOfDepVarsToSaveVector )
        throw std::runtime_error( "The specified save variable has not been initialized in the propagator settings" );

    positionInVector += indexInVectorialVariable_;

    return positionInVector;
}


double ObjectiveFunctionFromMinOrMaxDependentVariableSettings::getMinOrMaxSaveVariable(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator )
{
    int index = getPositionInDependentVariablesSaveVector( dynamicsSimulator );
    std::map< double, Eigen::VectorXd > dependentVariableHistory
            = dynamicsSimulator->getDependentVariableHistory();
    std::map< double, Eigen::VectorXd >::iterator it = dependentVariableHistory.begin();
    double bestValue = it->second[index];
    while( it != dependentVariableHistory.end() )
    {
        if( (minimumOrMaximum_ && ( it->second[index] < bestValue ))
                || (!minimumOrMaximum_ &&  (it->second[index] > bestValue )) )
            bestValue = it->second[index];
        ++it;
    }

    return bestValue;
}


}

}
