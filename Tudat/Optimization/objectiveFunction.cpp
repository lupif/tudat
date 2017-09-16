#include"objectiveFunction.h"


namespace tudat{

namespace optimization{

int ObjectiveFunctionFromFinalDependentVariableSettings::getPositionInDependentVariablesSaveVector(
        boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator ){

    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings >(
                dynamicsSimulator->getPropagatorSettings() );

    propagators::PropagationDependentVariables desiredDependentVariable =
            dependentVariableSaveSettings_->variableType_;

    if( indexInVectorialVariable_ >= propagators::getDependentVariableSize( desiredDependentVariable  ) )
        throw std::runtime_error( "The desired position in the chosed variable is higher than its size" );

    int sizeOfDepVarsToSaveVector = propagatorSettings->getDependentVariablesToSave()->dependentVariables_.size();
    int i = 0;
    int positionInVector = 0;
    bool found = false;
    do{
        propagators::PropagationDependentVariables currentDependentVariable =
                propagatorSettings->getDependentVariablesToSave()->dependentVariables_[i]->variableType_;
        if( currentDependentVariable == desiredDependentVariable )
        {
            if( desiredDependentVariable == propagators::single_acceleration_norm_dependent_variable ||
                    desiredDependentVariable == propagators::single_acceleration_dependent_variable )
            {
                boost::shared_ptr< propagators::SingleAccelerationDependentVariableSaveSettings > currentDependentVarSettings =
                        boost::dynamic_pointer_cast< propagators::SingleAccelerationDependentVariableSaveSettings >(
                            propagatorSettings->getDependentVariablesToSave()->dependentVariables_[i] );
                boost::shared_ptr< propagators::SingleAccelerationDependentVariableSaveSettings > desiredDependentVarSettings =
                            boost::dynamic_pointer_cast< propagators::SingleAccelerationDependentVariableSaveSettings >(
                                dependentVariableSaveSettings_);
                if( *desiredDependentVarSettings == *currentDependentVarSettings )
                    found = true;
            }
            else if( desiredDependentVariable == propagators::intermediate_aerodynamic_rotation_matrix_variable)
            {
                boost::shared_ptr< propagators::IntermediateAerodynamicRotationVariableSaveSettings > currentDependentVarSettings =
                        boost::dynamic_pointer_cast< propagators::IntermediateAerodynamicRotationVariableSaveSettings >(
                            propagatorSettings->getDependentVariablesToSave()->dependentVariables_[i] );
                boost::shared_ptr< propagators::IntermediateAerodynamicRotationVariableSaveSettings > desiredDependentVarSettings =
                            boost::dynamic_pointer_cast< propagators::IntermediateAerodynamicRotationVariableSaveSettings >(
                                dependentVariableSaveSettings_);
                if( *desiredDependentVarSettings == *currentDependentVarSettings )
                    found = true;
            }
            else if( desiredDependentVariable == propagators::relative_body_aerodynamic_orientation_angle_variable )
            {
                boost::shared_ptr< propagators::BodyAerodynamicAngleVariableSaveSettings > currentDependentVarSettings =
                    boost::dynamic_pointer_cast< propagators::BodyAerodynamicAngleVariableSaveSettings >(
                            propagatorSettings->getDependentVariablesToSave()->dependentVariables_[i] );
                boost::shared_ptr< propagators::BodyAerodynamicAngleVariableSaveSettings > desiredDependentVarSettings =
                    boost::dynamic_pointer_cast< propagators::BodyAerodynamicAngleVariableSaveSettings >(
                            dependentVariableSaveSettings_);
                if( *desiredDependentVarSettings == *currentDependentVarSettings )
                    found = true;
            }
            else if( *dependentVariableSaveSettings_ == propagatorSettings->getDependentVariablesToSave()->dependentVariables_[i] )
                found = true;
        }
        else
            positionInVector += propagators::propagators::getDependentVariableSize( currentDependentVariable );

        i++;
    }
    while( found == false || i < sizeOfDepVarsToSaveVector );

    if( found == false && i >= sizeOfDepVarsToSaveVector )
        throw std::runtime_error( "The specified save variable has not been initialized in the propagator settings" );

    positionInVector += indexInVectorialVariable_;

    return positionInVector;
}


double ObjectiveFunctionFromMinOrMaxDependentVariableSettings::getMinOrMaxSaveVariable(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator > dynamicsSimulator )
{
    int index = getPositionInDependentVariablesSaveVector( dynamicsSimulator );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator->getDependentVariableHistory();
    double bestValue = dependentVariableHistory.begin()->second[index];
    for( std::map< double, Eigen::VectorXd >::iterator it = dependentVariableHistory.begin();
            it!= dependentVariableHistory.end(); ++it ){
        if( (minimumOrMaximum_ && it->second[index] < bestValue)
                || (!minimumOrMaximum_ && it->second[index] > bestValue) )
            bestValue = it->second[index];
    }

    return bestValue;
}


}

}
