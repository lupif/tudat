/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include"missionSegmentSettings.h"

namespace tudat {

namespace optimization {


//Check whether the dynamics simulator has a certain type of propagator
bool MissionSegmentSettings::hasPropagatorSettings( propagators::IntegratedStateType propagatorType )
{

    bool found = false;

    boost::shared_ptr< propagators::PropagatorSettings< > > propagatorSettings =
            dynamicsSimulator_->getPropagatorSettings();

    if( propagatorSettings->stateType_ == propagatorType )
    {
        found = true;
    }
    else if( propagatorSettings->stateType_ == propagators::hybrid )
    {
        boost::shared_ptr< propagators::MultiTypePropagatorSettings< double > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< propagators::MultiTypePropagatorSettings< double > >(
                    propagatorSettings );
        std::map< propagators::IntegratedStateType,
                std::vector< boost::shared_ptr< propagators::PropagatorSettings< > > > > propagatorSettingsMap =
                multiTypePropagatorSettings->propagatorSettingsMap_;
        if( propagatorSettingsMap.count( propagatorType ) != 0 )
            found = true;
    }

    return found;
}


//Return the pointer to a certain type of propagator depending on the integrated state type
boost::shared_ptr< propagators::PropagatorSettings< > >  MissionSegmentSettings::selectPropagatorSettings(
        propagators::IntegratedStateType propagatorType )
{

    boost::shared_ptr< propagators::PropagatorSettings< > > propagatorSettings =
            dynamicsSimulator_->getPropagatorSettings();

    if( propagatorSettings->stateType_ == propagatorType )
    {
        return propagatorSettings;
    }
    else if( propagatorSettings->stateType_ == propagators::hybrid )
    {
        boost::shared_ptr< propagators::MultiTypePropagatorSettings< double > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< propagators::MultiTypePropagatorSettings< double > >(
                    propagatorSettings );
        std::map< propagators::IntegratedStateType,
                std::vector< boost::shared_ptr< propagators::PropagatorSettings< double > > > > propagatorSettingsMap =
                multiTypePropagatorSettings->propagatorSettingsMap_;
        if( propagatorSettingsMap.count( propagatorType ) != 0 )
            return propagatorSettingsMap[ propagatorType ][0];
        else
            throw std::runtime_error( "Mission linker: Unable to retrieve the propagator settings" );
    }
    else throw std::runtime_error( "Mission linker: Unable to retrieve the propagator settings" );
}


//Retrieve the initial orbital state associated to a body in the translational state propagator
//settings
Eigen::Vector6d MissionSegmentSettings::getInitialOrbitalState( const std::string selectedBodyName,
                                                                int* initialSegmentSize ){

    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            selectPropagatorSettings( propagators::transational_state ) );

    int index = -1;
    for( unsigned int i = 0; i<( propagatorSettings->bodiesToIntegrate_ ).size(); i++ )
    {
        if( selectedBodyName == propagatorSettings->bodiesToIntegrate_[i] ){
            index = i;
        }
        if( index == -1 )
          *initialSegmentSize = *initialSegmentSize + 6;
    }

    Eigen::VectorXd states = propagatorSettings->getInitialStates();
    Eigen::Vector6d actualState = states.segment( *initialSegmentSize, 6 );
    return actualState;

}


void MissionSegmentSettings::resetOrbitalState( Eigen::Vector6d modifiedState, int initialSegmentSize )
{

    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            selectPropagatorSettings( propagators::transational_state ) );
    Eigen::VectorXd bodyStates;
    bodyStates = propagatorSettings->getInitialStates();
    bodyStates.segment( initialSegmentSize, 6 ) = modifiedState;
    propagatorSettings->resetInitialStates( bodyStates );

    if( hasPropagatorSettings( propagators::hybrid ) )
    {
        dynamicsSimulator_->propagatorSettings_->resetInitialStates( getInitialStates() );
    }


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
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings  =
        boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings( propagators::transational_state ) );

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
        std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > >depVarSave =
                propagatorSettings->getDependentVariablesToSave()->dependentVariables_;

        double terminationValue = terminationSettings->limitValue_; //limit value
        int accumulator = 0;

        //Get the last entry in the history of dependent variables

        std::map< double, Eigen::VectorXd > dependentVariableHistory =
                dynamicsSimulator_->getDependentVariableHistory();

        std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                dependentVariableHistory.rbegin();

        for( unsigned int i = 0; i < depVarSave.size(); i++ )
        {
            // Look for a match between the dependent variable settings in the terminationSettings
            // and the dependent variable setting in the save variable list
           if( compareSingleDependentVariableSettings( terminationSettings->dependentVariableSettings_,
                                                       depVarSave[i]) )
           {

               // last time and value
               double time1 = rit->first;
               double value1 = rit->second[accumulator];
               ++rit;

               // penultimate time and value
               double time2 = rit->first;
               double value2 = rit->second[accumulator];

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
        std::map< double, Eigen::VectorXd > dependentVariableHistory =
                dynamicsSimulator_->getDependentVariableHistory();

        std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                dependentVariableHistory.rbegin();
        // scout all the termination settings in the hybrid termination conditions
        for( unsigned int i = 0; i <= hybridTerminationSettings->terminationSettings_.size(); i++ )

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
                if( rit->first >= timeTerminationSettings->terminationTime_ )
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
                std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > depVarSave =
                        propagatorSettings->getDependentVariablesToSave()->dependentVariables_;

                double terminationValue = terminationSettings->limitValue_; // limit value
                int accumulator = 0;

                // retrieve last entry in dependent variable history map
                std::map< double, Eigen::VectorXd > dependentVariableHistory =
                        dynamicsSimulator_->getDependentVariableHistory();

                std::map< double, Eigen::VectorXd >::reverse_iterator rit =
                        dependentVariableHistory.rbegin();

                // scout all the dependent variables to save
                for( unsigned int j = 0; j < depVarSave.size(); j++ )
                {
                   // look for a match in the dependent variable settings in the termination conditions
                   // and the dependend current variable to save
                   if( compareSingleDependentVariableSettings( terminationSettings->dependentVariableSettings_, depVarSave[j] ) )
                   {
                       // if the termination condition has been fulfilled
                       if( (terminationSettings->useAsLowerLimit_ && rit->second[accumulator] <= terminationSettings->limitValue_) ||
                            (!terminationSettings->useAsLowerLimit_ && rit->second[accumulator] >= terminationSettings->limitValue_) )
                       {
                            // final time and value
                            double time1 = rit->first;
                            double value1 = rit->second[accumulator];
                            ++rit;

                            // penultimate time and value
                            double time2 = rit->first;
                            double value2 = rit->second[accumulator];

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

std::map< std::string, Eigen::Vector6d > MissionSegmentSettings::getFinalOrbitalStates( const double finalTime )
{

    std::map< double, Eigen::VectorXd > numericalSolution =
            dynamicsSimulator_->getEquationsOfMotionNumericalSolution();
    std::map< double, Eigen::VectorXd >::reverse_iterator it =
            numericalSolution.rbegin();

    double time1 = it->first;
    Eigen::VectorXd states1 = it->second;
    ++it;
    double time2 = it->first;

    fflush(stdout);
    fflush(stdout);
    Eigen::VectorXd states2 = it->second;


    while( !( time1 >= finalTime && time2 <= finalTime ) ){

        time1 = time2;
        states1 = states2;
        ++it;
        time2 = it->first;
        states2 = it->second;

    }

    Eigen::VectorXd unsortedFinalStates = states2 + (states1 - states2)/(time1 - time2)*(time1 - finalTime);

    std::vector< std::string > bodies = boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings( propagators::transational_state ) )
            ->bodiesToIntegrate_;

    std::map< std::string, Eigen::Vector6d > finalStates;

    for( unsigned int i = 0; i < bodies.size(); i++ )
    {
        finalStates[ bodies[i] ] = unsortedFinalStates.segment( i*6, 6 );

    }

    return finalStates;


}

std::map< std::string, double > MissionSegmentSettings::getFinalBodyMasses( const double finalTime )
{

    std::map< double, Eigen::VectorXd > numericalSolution =
            dynamicsSimulator_->getEquationsOfMotionNumericalSolution();
    std::map< double, Eigen::VectorXd >::reverse_iterator it =
            numericalSolution.rbegin();

    double time1 = it->first;
    Eigen::VectorXd states1 = it->second;
    ++it;
    double time2 = it->first;
    Eigen::VectorXd states2 = it->second;

    while( !( time1 >= finalTime && time2 <= finalTime ) ){

        time1 = time2;
        states1 = states2;
        ++it;
        time2 = it->first;
        states2 = it->second;

    }

    Eigen::VectorXd unsortedFinalStates = states2 + (states1 - states2)/(time1 - time2)*(time1 - finalTime);

    std::vector< std::string > orbitalBodies = boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings( propagators::transational_state ) )
            ->bodiesToIntegrate_;

    const int initialMassesSegment = 6 * orbitalBodies.size();

    std::vector< std::string > bodies = boost::dynamic_pointer_cast< propagators::MassPropagatorSettings< double > >(
                selectPropagatorSettings( propagators::body_mass_state ) )
            ->bodiesWithMassToPropagate_;

    std::map< std::string, double > finalMasses;

    for( unsigned int i = 0; i < bodies.size(); i++ )
    {
        finalMasses[ bodies[i] ] = unsortedFinalStates( initialMassesSegment + i );
    }

    return finalMasses;

}



void MissionSegmentSettings::resetInitialOrbitalStates( std::map< std::string, Eigen::Vector6d > newStates ){

    std::vector< std::string > bodies = boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings( propagators::transational_state ) )->bodiesToIntegrate_;

    Eigen::VectorXd initialStates = boost::dynamic_pointer_cast<
            propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings( propagators::transational_state ) )->getInitialStates();

    for( unsigned int i = 0; i < bodies.size(); i++ )
    {
        initialStates.segment( i*6, 6 ) =   newStates[ bodies[i] ];
    }

    boost::dynamic_pointer_cast<
            propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings( propagators::transational_state ) )->resetInitialStates(initialStates);

}


void MissionSegmentSettings::resetStartEpoch( const double startEpoch ){

    dynamicsSimulator_->integratorSettings_->initialTime_ = startEpoch;

}


void MissionSegmentSettings::resetBodyMasses( std::map< std::string, double > newMasses )
{

    for( std::map< std::string, double >::iterator it = newMasses.begin(); it != newMasses.end(); ++it )
    {
        if( dynamicsSimulator_->bodyMap_.count( it->first ) > 0 )
            dynamicsSimulator_->bodyMap_[ it->first ]->setConstantBodyMass( it->second );
    }

    if( hasPropagatorSettings( propagators::body_mass_state ) )
    {

        boost::shared_ptr< propagators::MassPropagatorSettings< double > > propagatorSettings =
                boost::dynamic_pointer_cast< propagators::MassPropagatorSettings< double > >(
                    selectPropagatorSettings( propagators::body_mass_state ) );

        std::vector< std::string > bodies = propagatorSettings->bodiesWithMassToPropagate_;

        std::vector< double > masses;

        for( unsigned int i = 0; i < bodies.size(); i++ )
        {
            if( newMasses.count( bodies[i] ) > 0 )
                masses.push_back( newMasses[ bodies[i] ] );
        }

        Eigen::VectorXd newStates;
        newStates.resize( masses.size() );
        for( unsigned int i = 0; i < masses.size(); i++ )
            newStates[i] = masses[i];

        propagatorSettings->resetInitialStates( newStates );

    }

}



void MissionSegmentSettings::setPreliminaryConditions( const double startEpoch,
                             std::map< std::string, Eigen::Vector6d > newStates,
                             std::map< std::string, double > newMasses ){

    resetInitialOrbitalStates( newStates );
    resetBodyMasses( newMasses );

    if( hasPropagatorSettings( propagators::hybrid ) )
    {
        boost::dynamic_pointer_cast< propagators::MultiTypePropagatorSettings< double > >
                ( selectPropagatorSettings( propagators::hybrid ) )->resetInitialStates(
                    getInitialStates() );
    }
    resetStartEpoch( startEpoch );


}


void MissionSegmentSettings::setPropagationConditions( Eigen::VectorXd& newValues )
{
    int accumulator = 0;
    for( unsigned int i = 0; i < decisionVariableSettings_->decisionVariableSettings_.size(); i++ ){

        if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ == simulation_time_decision_variable )
        {
           boost::shared_ptr< propagators::PropagatorSettings< > > propagatorSettings =
                dynamicsSimulator_->getPropagatorSettings();

           boost::shared_ptr< propagators::PropagationTimeTerminationSettings > timeTerminationSettings;
           if( propagatorSettings->getTerminationSettings()->terminationType_ ==
                   propagators::time_stopping_condition )
               timeTerminationSettings = boost::dynamic_pointer_cast<
                       propagators::PropagationTimeTerminationSettings >( propagatorSettings->getTerminationSettings());
           else if( propagatorSettings->getTerminationSettings()->terminationType_ ==
                    propagators::hybrid_stopping_condition )
           {
               boost::shared_ptr< propagators::PropagationHybridTerminationSettings > hybridTerminationSettings =
                       boost::dynamic_pointer_cast< propagators::PropagationHybridTerminationSettings >(
                           propagatorSettings->getTerminationSettings() );
               for( unsigned int j = 0; j < hybridTerminationSettings->terminationSettings_.size(); j++)
               {
                   if(hybridTerminationSettings->terminationSettings_[j]->terminationType_
                           == propagators::time_stopping_condition)
                       timeTerminationSettings = boost::dynamic_pointer_cast<
                               propagators::PropagationTimeTerminationSettings >(
                                       hybridTerminationSettings->terminationSettings_[j] );
               }

           }
           if( timeTerminationSettings == NULL )
               throw std::runtime_error( "Error: no time termination settings found" );
           else
               timeTerminationSettings->terminationTime_ =
                       dynamicsSimulator_->getIntegratorSettings()->initialTime_ + newValues[accumulator];
           accumulator++;
        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_
                == from_termination_settings_decision_variable )
        {
            boost::shared_ptr< SingleDecisionVariableFromTerminationSettings > decisionVariableSettings =
                   boost::dynamic_pointer_cast< SingleDecisionVariableFromTerminationSettings >(
                       decisionVariableSettings_->decisionVariableSettings_[i]);
            *(decisionVariableSettings->memoryPositionOfVariable_) = newValues[accumulator];
            accumulator++;

        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ ==
                initial_cartesian_state_decision_variable )
        {
            int initialSegmentSize = 0;
            getInitialOrbitalState( decisionVariableSettings_->decisionVariableSettings_[i]->associatedBody_,
                                    &initialSegmentSize );
            resetOrbitalState(newValues.segment(accumulator, 6), initialSegmentSize );
            accumulator += 6;
        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ ==
                initial_cartesian_velocity_decision_variable )
        {
            int initialSegmentSize = 0;
            Eigen::Vector6d cartesianState =  getInitialOrbitalState(
                        decisionVariableSettings_->decisionVariableSettings_[i]->associatedBody_,
                                &initialSegmentSize );
            cartesianState.segment(3, 3) = newValues.segment(accumulator, 3);
            resetOrbitalState( cartesianState, initialSegmentSize );
            accumulator += 3;
        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ ==
                initial_cartesian_position_decision_variable )
        {
            int initialSegmentSize = 0;
            Eigen::Vector6d cartesianState =  getInitialOrbitalState(
                        decisionVariableSettings_->decisionVariableSettings_[i]->associatedBody_,
                                &initialSegmentSize );
            cartesianState.segment(0, 3) = newValues.segment(accumulator, 3);
            resetOrbitalState(cartesianState, initialSegmentSize );
            accumulator += 3;
        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ ==
                 single_cartesian_component_decision_variable )
        {
            boost::shared_ptr< SingleCartesianComponentDecisionVariableSettings > decisionVariableSettings
                    = boost::dynamic_pointer_cast< SingleCartesianComponentDecisionVariableSettings >(
                        decisionVariableSettings_->decisionVariableSettings_[i]);
            modifyCartesianComponent( decisionVariableSettings->cartesianComponent_, decisionVariableSettings->associatedBody_,
                                      newValues[accumulator] );
            accumulator++;
        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ ==
                 single_kepler_element_decision_variable )
        {
            boost::shared_ptr< SingleKeplerElementDecisionVariableSettings > decisionVariableSettings
                    = boost::dynamic_pointer_cast< SingleKeplerElementDecisionVariableSettings >(
                        decisionVariableSettings_->decisionVariableSettings_[i]);
            modifyKeplerElement( decisionVariableSettings->keplerElement_, decisionVariableSettings->associatedBody_,
                                 decisionVariableSettings->centralBody_, newValues[accumulator] );
            accumulator++;
        }
        else if( decisionVariableSettings_->decisionVariableSettings_[i]->decisionVariable_ ==
                 single_spherical_orbital_element_decision_variable )
        {
            boost::shared_ptr< SingleSphericalOrbitalElementDecisionVariableSettings > decisionVariableSettings
                    = boost::dynamic_pointer_cast< SingleSphericalOrbitalElementDecisionVariableSettings >(
                        decisionVariableSettings_->decisionVariableSettings_[i]);
            modifySphericalComponent( decisionVariableSettings->sphericalOrbitalElement_, decisionVariableSettings->associatedBody_,
                                      newValues[accumulator] );
            accumulator++;
        }

    }

    dynamicsSimulator_->resetPropagationTerminationConditions();

}


//Return the value of the objective function with the current configuration
//This function is to be run after the propagation has taken place.
double MissionSegmentSettings::getObjectiveFunctionValue( void )
{

    // Get the final epoch interpolated according to the termination settings.
    double finalTime = getFinalSimulationEpoch();

    double objectiveFunctionValue;

    //Retrieve a final cartesian component
    if( objectiveFunctionSettings_->objectiveVariable_ == final_cartesian_component_objective_function )
    {
        // Dynamic cast for objective function settings
        boost::shared_ptr< FinalCartesianComponentObjectiveFunctionSettings > objectiveFunctionSettings =
                boost::dynamic_pointer_cast< FinalCartesianComponentObjectiveFunctionSettings >(
                    objectiveFunctionSettings_);
        // Get the final orbital state interpolated to the final epoch
        std::map< std::string, Eigen::Vector6d > finalStates = getFinalOrbitalStates( finalTime );
        objectiveFunctionValue =  finalStates[objectiveFunctionSettings->associatedBody_][
                objectiveFunctionSettings->cartesianComponent_];
    }

    else if( objectiveFunctionSettings_->objectiveVariable_ == final_kepler_orbital_element_objective_function )
    {
        boost::shared_ptr< FinalKeplerElementObjectiveFunctionSettings > objectiveFunctionSettings =
                boost::dynamic_pointer_cast< FinalKeplerElementObjectiveFunctionSettings >(
                    objectiveFunctionSettings_);
        std::map< std::string, Eigen::Vector6d > finalStates = getFinalOrbitalStates( finalTime );
        Eigen::Vector6d keplerState = orbital_element_conversions::convertCartesianToKeplerianElements(
                    finalStates[objectiveFunctionSettings->associatedBody_], dynamicsSimulator_->bodyMap_[
                objectiveFunctionSettings->centralBody_]->getGravityFieldModel()->getGravitationalParameter());
        objectiveFunctionValue = keplerState[objectiveFunctionSettings->keplerElement_];
    }
    else if( objectiveFunctionSettings_->objectiveVariable_ == final_spherical_orbital_element_objective_function )
    {
        boost::shared_ptr< FinalSphericalOrbitalElementObjectiveFunctionSettings > objectiveFunctionSettings =
                boost::dynamic_pointer_cast< FinalSphericalOrbitalElementObjectiveFunctionSettings >(
                    objectiveFunctionSettings_);
        std::map< std::string, Eigen::Vector6d > finalStates = getFinalOrbitalStates( finalTime );
        Eigen::Vector6d sphericalState = orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                    finalStates[objectiveFunctionSettings->associatedBody_] );
        objectiveFunctionValue = sphericalState[objectiveFunctionSettings->sphericalOrbitalElement_];
    }
    else if( objectiveFunctionSettings_->objectiveVariable_ == final_dependent_variable_objective_function )
    {
        boost::shared_ptr< ObjectiveFunctionFromFinalDependentVariableSettings > objectiveFunctionSettings =
                boost::dynamic_pointer_cast< ObjectiveFunctionFromFinalDependentVariableSettings >(
                    objectiveFunctionSettings_ );
        int positionInVector = objectiveFunctionSettings->getPositionInDependentVariablesSaveVector(
                    dynamicsSimulator_ );
        std::map< double, Eigen::VectorXd > dependentVariableHistory
                = dynamicsSimulator_->getDependentVariableHistory();
        std::map< double, Eigen::VectorXd >::reverse_iterator it = dependentVariableHistory.rbegin();

        double time1 = it->first;
        double value1 = it->second[positionInVector];

        it++;

        double time2 = it->first;
        double value2 = it->second[positionInVector];

        while( !( time2 <= finalTime && time1 >= finalTime) )
        {
            time1 = time2;
            value1 = value2;

            it++;

            time2 = it->first;
            value2 = it->second[positionInVector];

        }

        objectiveFunctionValue = value2 + ( value1 - value2 )*(finalTime - time2)/(time1 - time2);
    }
    else if( objectiveFunctionSettings_->objectiveVariable_ == min_max_dependent_variable_objective_function )
    {
        boost::shared_ptr< ObjectiveFunctionFromMinOrMaxDependentVariableSettings > objectiveFunctionSettings =
                boost::dynamic_pointer_cast< ObjectiveFunctionFromMinOrMaxDependentVariableSettings >(
                    objectiveFunctionSettings_ );

        objectiveFunctionValue = objectiveFunctionSettings->getMinOrMaxSaveVariable(dynamicsSimulator_);

    }
    else if( objectiveFunctionSettings_->objectiveVariable_ == user_defined_objective_function )
    {
        boost::shared_ptr< UserDefinedObjectiveFunctionSettings > objectiveFunctionSettings =
                boost::dynamic_pointer_cast< UserDefinedObjectiveFunctionSettings >(
                    objectiveFunctionSettings_ );
        objectiveFunctionValue = objectiveFunctionSettings->userDefinedFunction_();
    }

    if( objectiveFunctionSettings_->objectiveValueIsSet_ == true )
        return fabs( objectiveFunctionSettings_->objectiveValue_ - objectiveFunctionValue );
    else if( objectiveFunctionSettings_->minimizeOrMaximize_ )
        return objectiveFunctionValue;
    else
        return -objectiveFunctionValue;

    return -1.0;

}


Eigen::VectorXd MissionSegmentSettings::getInitialStates( void )
{

    Eigen::VectorXd initialStates;

    initialStates = boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
                selectPropagatorSettings(propagators::transational_state ) )->getInitialStates();

    if( hasPropagatorSettings( propagators::body_mass_state ) )
    {
        boost::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings =
                boost::dynamic_pointer_cast< propagators::MassPropagatorSettings< double > >(
                                selectPropagatorSettings( propagators::body_mass_state ) );

        Eigen::VectorXd masses = massPropagatorSettings->getInitialStates();

        Eigen::VectorXd newInitialStates = Eigen::VectorXd::Zero( initialStates.rows() + masses.rows() );

        newInitialStates.segment(0, initialStates.rows() ) = initialStates;

        newInitialStates.segment( initialStates.rows(), masses.rows() ) = masses;

        return newInitialStates;

    }

    return initialStates;


}



}


}
