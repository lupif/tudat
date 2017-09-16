#include"missionSegmentSettings.h"


namespace tudat {

namespace optimization {


Eigen::Vector6d MissionSegmentSettings::getInitialOrbitalState( const std::string selectedBodyName,
                                                                int* initialSegmentSize ){
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            dynamicsSimulator_->getPropagatorSettings() );
    int index = -1;
    for( unsigned int i = 0; i<( propagatorSettings->bodiesToIntegrate_ ).size(); i++ )
    {
        if( selectedBodyName == propagatorSettings->bodiesToIntegrate_[i] ){
            index = i;
        }
        if( index == -1 )
          *initialSegmentSize = *initialSegmentSize + 6;
    }

    Eigen::VectorXd states;
    Eigen::Vector6d actualState;
    states.resize( propagatorSettings->getInitialStates().rows() );
    states = propagatorSettings->getInitialStates();
    actualState = states.segment( *initialSegmentSize, 6 );
    return actualState;
}


void MissionSegmentSettings::resetOrbitalState( Eigen::Vector6d& modifiedState, int initialSegmentSize )
{
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            dynamicsSimulator_->getPropagatorSettings() );
    Eigen::VectorXd bodyStates;
    bodyStates = propagatorSettings->getInitialStates();
    bodyStates.segment( initialSegmentSize, 6 ) = modifiedState;
    propagatorSettings->resetInitialStates( bodyStates );

}

void MissionSegmentSettings::modifyCartesianComponent( orbital_elements::CartesianElements cartesianComponent,
                               const std::string selectedBodyName, const double newValue )
{
      int initialSegmentSize = 0;
      Eigen::Vector6d cartesianState = getInitialOrbitalState(selectedBodyName,
                                                              &initialSegmentSize);
      cartesianState[ cartesianComponent ] = newValue;
      resetOrbitalState( cartesianState, initialSegmentSize );
}


void MissionSegmentSettings::modifyKeplerElement(orbital_elements::KeplerianElements keplerElement,
                         const std::string selectedBodyName, const std::string centralBody, const double newValue )
{

    int initialSegmentSize = 0;
    Eigen::Vector6d cartesianState = getInitialOrbitalState(selectedBodyName,
                                                            &initialSegmentSize);
    Eigen::Vector6d keplerState = orbital_element_conversions::convertCartesianToKeplerianElements( cartesianState,
            dynamicsSimulator_->getNamedBodyMap()[centralBody]->getGravityFieldModel()->getGravitationalParameter() );
    keplerState[keplerElement] = newValue;
    cartesianState =  orbital_element_conversions::convertKeplerianToCartesianElements( keplerState,
            dynamicsSimulator_->getNamedBodyMap()[centralBody]->getGravityFieldModel()->getGravitationalParameter());
    resetOrbitalState( cartesianState, initialSegmentSize );

}


void MissionSegmentSettings::modifySphericalComponent( orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
                               const std::string selectedBodyName, const double newValue )
{
    int initialSegmentSize = 0;
    Eigen::Vector6d cartesianState = getInitialOrbitalState( selectedBodyName,
                                                            &initialSegmentSize );
    Eigen::Vector6d sphericalState = orbital_element_conversions::convertCartesianToSphericalOrbitalState( cartesianState );
    sphericalState[ sphericalOrbitalElement ] = newValue;
    cartesianState = orbital_element_conversions::convertSphericalOrbitalToCartesianState( sphericalState );
    resetOrbitalState( cartesianState, initialSegmentSize );

}


double MissionSegmentSettings::getFinalSimulationEpoch(){

    // Retrieve the translationalStatePropagatorSettings from the dynamics simulator
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings > propagatorSettings  =
        boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            dynamicsSimulator_->getPropagatorSettings() );

    // Get the termination settings type
    propagators::PropagationTerminationTypes terminationType = propagatorSettings->getTerminationSettings()->
            terminationType_;

    // If time stopping then return the conditions for time stopping
    if( terminationType == propagators::time_stopping_condition )
        return boost::dynamic_pointer_cast< propagators::PropagationTimeTerminationSettings >(
                    propagatorSettings->getTerminationSettings() )->terminationTime_;
    // If dependent variable stopping condition
    else if( terminationType == propagators::dependent_variable_stopping_condition )
    {
        // Retrieve the variable termination settings
        boost::shared_ptr< propagators::PropagationDependentVariableTerminationSettings > terminationSettings =
                boost::dynamic_pointer_cast< propagators::PropagationDependentVariableTerminationSettings >(
                    propagatorSettings->getTerminationSettings() );
        // Retrieve the saved variables
        std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > depVarSave =
                propagatorSettings->getDependentVariablesToSave()->dependentVariables_;

        double terminationValue = terminationSettings->limitValue_; //limit value
        int accumulator = 0;

        //Get the last entry in the history of dependent variables
        std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                dynamicsSimulator_->getDependentVariableHistory().rbegin();

        for( int i = 0; i < depVarSave.size(); i++ )
        {
            // Look for a match between the dependent variable settings in the terminationSettings
            // and the dependent variable setting in the save variable list
           if( *(terminationSettings->dependentVariableSettings_) == *(depVarSave[i]) )
           {

               // last time and value
               double time1 = *rit->first;
               double value1 = *rit->second[accumulator];
               ++rit;

               // penultimate time and value
               double time2 = *rit->first;
               double value2 = *rit->second[accumulator];

               // interpolate time
               double finalTime = time2 + ( time1 - time2 )/(value1 - value2)*(terminationValue - value2);
               return finalTime;
           }
           else
               // go to the next dependent variable
               accumulator += propagators::getDependentVariableSize( depVarSave[i]->variableType_ );
        }
    }
    // if hybrid conditions:
    else if( terminationType == propagators::hybrid_stopping_condition )
    {
        // retrieve hybrid termination conditions
        boost::shared_ptr< propagators::PropagationHybridTerminationSettings > hybridTerminationSettings =
                boost::dynamic_pointer_cast< propagators::PropagationHybridTerminationSettings >(
                    propagatorSettings->getTerminationSettings() );
        // retrieve last entry in map of dependent variable history
        std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                dynamicsSimulator_->getDependentVariableHistory().rbegin();
        // scout all the termination settings in the hybrid termination conditions
        for( int i = 0; i <= hybridTerminationSettings->terminationSettings_.size(); i++ )

            // if termination time settings
            if( hybridTerminationSettings->terminationSettings_[i]->terminationType_ ==
                    propagators::time_stopping_condition )
            {
                // retrieve time termination settings
                boost::shared_ptr< propagators::PropagationTimeTerminationSettings > timeTerminationSettings =
                        boost::dynamic_pointer_cast< propagators::PropagationTimeTerminationSettings >(
                                hybridTerminationSettings->terminationSettings_[i] );
                // if condition has been fulfilled, simply return the desired time termination
                // limit value
                if( *rit->first >= timeTerminationSettings->terminationTime_ )
                    return timeTerminationSettings->terminationTime_;
            }
            // if dependent variable termination settings
            else if( hybridTerminationSettings->terminationSettings_[i]->terminationType_ ==
                     propagators::dependent_variable_stopping_condition )
            {
                // retrieve dependent variable termination settings
                boost::shared_ptr< propagators::PropagationDependentVariableTerminationSettings > terminationSettings =
                        boost::dynamic_pointer_cast< propagators::PropagationDependentVariableTerminationSettings >(
                            propagatorSettings->getTerminationSettings() );
                // retrieve map of dependent variable history
                std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > depVarSave =
                        propagatorSettings->getDependentVariablesToSave()->dependentVariables_;

                double terminationValue = terminationSettings->limitValue_; // limit value
                int accumulator = 0;

                // retrieve last entry in dependent variable history map
                std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                        dynamicsSimulator_->getDependentVariableHistory().rbegin();

                // scout all the dependent variables to save
                for( int j = 0; j < depVarSave.size(); j++ )
                {
                   // look for a match in the dependent variable settings in the termination conditions
                   // and the dependend current variable to save
                   if( *(terminationSettings->dependentVariableSettings_) == *(depVarSave[j]) )
                   {
                       // if the termination condition has been fulfilled
                       if( (terminationSettings->useAsLowerLimit_ && *rit->second[accumulator] <= terminationSettings->limitValue_) ||
                            (!terminationSettings->useAsLowerLimit_ && *rit->second[accumulator] >= terminationSettings->limitValue_) )
                       {
                            // final time and value
                            double time1 = *rit->first;
                            double value1 = *rit->second[accumulator];
                            ++rit;

                            // penultimate time and value
                            double time2 = *rit->first;
                            double value2 = *rit->second[accumulator];

                            // interpolate time to match the desired termination value and return
                            double finalTime = time2 + ( time1 - time2 )/(value1 - value2)*(terminationValue - value2);
                            return finalTime;
                       }
                   }
                   else
                       // otherwise progress in the vector of saved variables
                       accumulator += propagators::getDependentVariableSize( depVarSave[j]->variableType_ );
                }
            }
        }
    // standard return value
    return -1.0;
}

}


}
