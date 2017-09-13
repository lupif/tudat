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
    case simulationTime:
        size = 1;
        break;
    case initialStateCartesianComponents:
        size = 6;
        break;
    case initialVelocitytateCartesianComponents:
        size = 3;
        break;
    case initialPositionStateCartesianComponents:
        size = 3;
        break;
    case singleCartesianComponentsElement:
        size = 1;
        break;
    case singleKeplerOrbitalElement:
        size = 1;
        break;
    case singleSphericalOrbitalElement:
        size = 1;
        break;
    case fromTerminationSettings:
        break;
    default:
        std::string errorMessage = "Error, did not recognize decision variable size of type: " +
                boost::lexical_cast< std::string >( decisionVariable );
        throw std::runtime_error( errorMessage );
    }

    return size;
}



SingleDecisionVariableSettings::SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                boost::shared_ptr< Boundaries > boundaries ) :
    decisionVariable_(decisionVariable), boundaries_( boundaries ){ }


SingleDecisionVariableSettings::SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                double lowerBoundary, double upperBoundary ) :
    decisionVariable_(decisionVariable)
{
    if( getDecisionVariableSize( decisionVariable ) != 1 && decisionVariable != fromTerminationSettings ){
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
                                Eigen::VectorXd& lowerBoundary, Eigen::VectorXd& upperBoundary) :
    decisionVariable_(decisionVariable)
{
    if( lowerBoundary.rows() != upperBoundary.rows() ){
        throw std::runtime_error( "Lower boundary and upper boundary size mismatch" );
    }
    else if( getDecisionVariableSize( decisionVariable ) != lowerBoundary.rows()
             && decisionVariable != fromTerminationSettings ){
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


Eigen::Vector6d getInitialOrbitalState( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                                        const std::string& selectedBodyName, int* initialSegmentSize ){
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            dynamicsSimulator->getPropagatorSettings() );
    int index = -1;
    int segmentSize = 0;
    for( unsigned int i = 0; i<( propagatorSettings->bodiesToIntegrate_ ).size(); i++ )
    {
        if( selectedBodyName == propagatorSettings->bodiesToIntegrate_[i] ){
            index = i;
            segmentSize = propagatorSettings->getStateSize();
        }
        if( index == -1 )
          *initialSegmentSize = *initialSegmentSize + propagatorSettings->getStateSize();
    }

    if( segmentSize != 6 )
        throw std::runtime_error("The defined state to modify is not an orbital state (6 components)");

    Eigen::VectorXd states;
    Eigen::VectorXd actualState;
    states.resize( propagatorSettings->getInitialStates().rows() );
    states = propagatorSettings->getInitialStates();
    actualState.resize( segmentSize );
    actualState = states.segment( *initialSegmentSize, segmentSize );
    return actualState;
}


void resetOrbitalState( boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                        Eigen::Vector6d& modifiedState, int initialSegmentSize )
{
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings =
            boost::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< > >(
            dynamicsSimulator->getPropagatorSettings() );
    Eigen::VectorXd bodyStates;
    bodyStates = propagatorSettings->getInitialStates();
    bodyStates.segment( initialSegmentSize, 6 ) = modifiedState;
    propagatorSettings->resetInitialStates( bodyStates );

}

void modifyCartesianComponent( orbital_elements::CartesianElements cartesianComponent,
                               boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                               const std::string& selectedBodyName, const double newValue )
{
      int initialSegmentSize = 0;
      Eigen::Vector6d cartesianState = getInitialOrbitalState(dynamicsSimulator, selectedBodyName,
                                                              &initialSegmentSize);
      cartesianState[ cartesianComponent ] = newValue;
      resetOrbitalState(dynamicsSimulator, cartesianState, initialSegmentSize );
}


void modifyKeplerElement(orbital_elements::KeplerianElements keplerElement,
                         boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                         const std::string& selectedBodyName, const std::string& centralBody, const double newValue )
{

    int initialSegmentSize = 0;
    Eigen::Vector6d cartesianState = getInitialOrbitalState(dynamicsSimulator, selectedBodyName,
                                                            &initialSegmentSize);
    Eigen::Vector6d keplerState = orbital_element_conversions::convertCartesianToKeplerianElements( cartesianState,
            dynamicsSimulator->getNamedBodyMap()[centralBody]->getGravityFieldModel()->getGravitationalParameter() );
    keplerState[keplerElement] = newValue;
    cartesianState =  orbital_element_conversions::convertKeplerianToCartesianElements( keplerState,
            dynamicsSimulator->getNamedBodyMap()[centralBody]->getGravityFieldModel()->getGravitationalParameter());
    resetOrbitalState(dynamicsSimulator, cartesianState, initialSegmentSize );

}


void modifySphericalComponent( orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
                               boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
                               const std::string& selectedBodyName, const double newValue )
{
    int initialSegmentSize = 0;
    Eigen::Vector6d cartesianState = getInitialOrbitalState( dynamicsSimulator, selectedBodyName,
                                                            &initialSegmentSize );
    Eigen::Vector6d sphericalState = orbital_element_conversions::convertCartesianToSphericalOrbitalState( cartesianState );
    sphericalState[ sphericalOrbitalElement ] = newValue;
    cartesianState = orbital_element_conversions::convertSphericalOrbitalToCartesianState( sphericalState );
    resetOrbitalState( dynamicsSimulator, cartesianState, initialSegmentSize );

}


} //

}
