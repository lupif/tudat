/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include"missionLinker.h"

#include"pagmo/algorithms/de1220.hpp"
#include"pagmo/algorithms/compass_search.hpp"
#include"pagmo/island.hpp"
#include"pagmo/problem.hpp"


namespace tudat{

namespace optimization{


//This function defines how many scalar decision variables are in each segment
//by filling a vector containing the same number of elements as the number of
//missionSegmentSettings with the integer number
/*!
 * Retrieve the number of decision variables for each mission segment
 * \param missionSegments Vector of mission segments
 * \return Vector of the same size as the input vector with the associated number
 * of decision variables for each mission segment
 */
std::vector< int > getNumberOfDecisionVariables(
        std::vector< boost::shared_ptr< MissionSegmentSettings > >& missionSegments )
{

    std::vector< int > nOfDecisionVariables; //return variable
    //Scout all the missionSegments
    for( unsigned int k = 0; k < missionSegments.size(); k++ )
    {
        nOfDecisionVariables.push_back(0); // standard is 0
        //If the decisionVariableSettings_ member is set, scout all the
        //singleDecisionVariableSetting in the vector inside
        //decisionVariableSettings and simply sum all the rows in the lower boundary
        if( missionSegments[k]->decisionVariableSettings_ != NULL )
        {
            std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > singleDecisionVarSettingsVec =
                    missionSegments[k]->decisionVariableSettings_->decisionVariableSettings_;
            for( unsigned int i = 0; i<singleDecisionVarSettingsVec.size(); i++ )
            {
                for( unsigned int j = 0; j < singleDecisionVarSettingsVec[i]->boundaries_->
                        getLowerBoundary().rows(); j++ ){
                    nOfDecisionVariables.back() += 1;
                }
            }
        }
    }
    return nOfDecisionVariables;

}

//This function retrieves all the boundaries in the missionSegments decision variables
//And stores them in a pair of vectors, as for Pagmo requirements.
/*!
 * Returns a pair of vectors containing respectively the lower and upper boundaries of all the decision variable
 * settings contained in a vector of MissionSegmentSetting objects
 * \param missionSegments Vector of MissionSegmentSetting objects
 * \return pair of vectors containing respectively the lower and upper boundaries of all the decision variable
 * settings contained in the parameter missionSegments
 */
std::pair< std::vector< double >, std::vector< double > > getBoundariesPair(
        std::vector< boost::shared_ptr< MissionSegmentSettings > >& missionSegments )
{
    std::pair< std::vector< double >, std::vector< double > > boundaries; //return variable

    //Scout all the mission segments
    for( unsigned int k = 0; k < missionSegments.size(); k++ )
    {
        //If decisionVariableSettings_ is set
        if( missionSegments[k]->decisionVariableSettings_ != NULL )
        {
            std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > singleDecisionVarSettingsVec =
                    missionSegments[k]->decisionVariableSettings_->decisionVariableSettings_;
            //Scout all the singleDecisionVariableSettings in the vector
            for( unsigned int i = 0; i<singleDecisionVarSettingsVec.size(); i++ )
            {
                //Transfer the values.
                for( unsigned int j = 0; j < singleDecisionVarSettingsVec[i]->boundaries_->
                        getLowerBoundary().rows(); j++ ){
                    boundaries.first.push_back(
                                singleDecisionVarSettingsVec[i]->boundaries_->getLowerBoundary()[j] );
                    boundaries.second.push_back(
                                singleDecisionVarSettingsVec[i]->boundaries_->getUpperBoundary()[j] );
                }
            }
        }
    }
    return boundaries;

}

//This function recovers all the singleDecisionVariableSettings objects in all the mission
//segments and stores them in a vector
//! Get all the single decision variable settings from a vector of mission segments
/*!
 * Retrieves a vector containing all the SingleDecisionVariableSettings objects from a vector of MissionSegmentSettings
 * objects.
 * \param missionSegmentsBlock Vector fo MissionSegmentSettings objects
 * \return A vector containing all the SingleDecisionVariableSettings objects from missionSegmentsBlock
 */
std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > getAllSingleDecisionVariableSettings(
        std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentsBlock )
{

    std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > decisionVariableSettings;
    for( unsigned int i = 0; i < missionSegmentsBlock.size(); i++)
    {

        if( missionSegmentsBlock[i]->decisionVariableSettings_ != NULL)
            for( unsigned int j = 0; j < missionSegmentsBlock[i]->decisionVariableSettings_
                    ->decisionVariableSettings_.size(); j++ )
                decisionVariableSettings.push_back( missionSegmentsBlock[i]->decisionVariableSettings_
                        ->decisionVariableSettings_[j] );
    }

    return decisionVariableSettings;
}

//This function retrievesl all the boundaries in a vector of SingleDecisionVariableSettings
//and stores them in a vector of pairs od Eigen::VectorXd. Each pair belongs to a
//singleDecisionVariableSettigns
//! Retrieve all the boundaries from a vector of SingleDecisionVariableSettings
/*!
 * Retrieve all the boundaries from a vector of SingleDecisionVariableSettings. Each multivariate boundary
 * corresponding to a SingleDecisionVariableSettings object in the vector is stored in a pair of Eigen.VectorXd.
 * The function returns a vector of pairs of the same size as the input vector.
 * \param decisionVariableSettings Vector of SingleDecisionVariableSettings objects
 * \return Vector of pairs of Eigen.VectorXd containing respectively the multivariate lower and upper boundaries.
 * Each pair corresponds to a DecisionVariableSettings object in decisionVariableSettings.
 */
std::vector< std::pair< Eigen::VectorXd, Eigen::VectorXd > > getAllBoundaries(
        std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > decisionVariableSettings )
{

    std::vector< std::pair< Eigen::VectorXd, Eigen::VectorXd > > boundariesVector;
    for( unsigned int i = 0; i < decisionVariableSettings.size(); i++ )
    {
        std::pair< Eigen::VectorXd, Eigen::VectorXd > boundaries;
        boundaries.first = decisionVariableSettings[i]->boundaries_->getLowerBoundary();
        boundaries.second = decisionVariableSettings[i]->boundaries_->getUpperBoundary();
        boundariesVector.push_back( boundaries );
    }

    return boundariesVector;

}

//! Set the multivariate boundaries in a vector of SingleDecisionVariableSettings objects
/*!
 * Set the multivariate boundaries of a vector of SingleDecisionVariableSettings obejcts.
 * \param boundariesVector Vector of pairs of Eigen.VectorXd being respectively the new values for the
 * lower and upper boundaries.
 * \param decisionVariableSettings vector of SingleDecisionVariableSettings objects whose boundaries
 * need to be reset.
 */
void resetBoundaries(  std::vector< std::pair< Eigen::VectorXd, Eigen::VectorXd > > boundariesVector,
                       std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > decisionVariableSettings )
{

    for( unsigned int i = 0; i < boundariesVector.size(); i++ )
    {
        decisionVariableSettings[i]->boundaries_->setLowerBoundary( boundariesVector[i].first );
        decisionVariableSettings[i]->boundaries_->setLowerBoundary( boundariesVector[i].second );
    }
}


// Pagmo problem used ot link the mission segments.
// The purpose of this problem is to retrieve the objective function value at the
// end of a chain of mission segments. The chain is defined as a series of mission
// segments with only one objective function.
//! Pagmo problem used ot link the mission segments.
struct LinkingProblem{

    // Empty constructor.
    // In its absence multithreading fail to recognize LinkingProblem as a Pagmo Problem type
    LinkingProblem( ){ }

    // Constructor defines the decision variables  boundary vectors and the segmentation associated
    // to each mission segment
    //! Constructor
    /*!
     * Constructor of the Pagmo problem. It accepts a set of mission segments and stores useful information
     * for the optimization, namely the number of single decision variables for each segment and
     * the Pagmo boundaries from all the missionSegments.
     * \param missionSegments Vector of MissionSegmentSettings whose first object needs to have
     * a DecisionVariableSettings object and last object needs to have an ObjectiveFunctionSettings object.
     */
    LinkingProblem( std::vector< boost::shared_ptr< MissionSegmentSettings > >& missionSegments ) :
        missionSegments_( missionSegments )
    {
        nOfDecisionVariables_ = getNumberOfDecisionVariables( missionSegments_ );
        boundaries_ = getBoundariesPair( missionSegments_ );
    }

    ~LinkingProblem( ){ }

    // Calculate the fitness, i.e. the objective function of the last mission segment
    // in the chain
    //! Calculate the fitness of the problem as per Pagmo standard
    /*!
     * Calculate the fitness, i.e. the value of the objective function at the end of the block
     * of mission segments.
     * \param decisionVariables Values of the decision variables of all the missionSegments in the vector
     * missionSegments_
     * \return Value of the objective function
     */
    std::vector< double > fitness( const std::vector< double > &decisionVariables ) const
    {
        int accumulator = 0;
        Eigen::VectorXd newValues;

        for( unsigned int i = 0; i < missionSegments_.size(); i++ ){

            // Set the (interpolated) final states of the previous mission segment
            if( i > 0 ){

                // get the interpolated end of simulation's epoch
                double prevFinalTime = missionSegments_[i-1]->getFinalSimulationEpoch();

                // get the interpolated final state of the simulation's epoch
                std::map< std::string, Eigen::Vector6d > prevFinalStates = missionSegments_[i-1]->
                        getFinalOrbitalStates( prevFinalTime );

                std::map< std::string, double > prevFinalMasses = std::map< std::string, double >();


                if( missionSegments_[i-1]->hasPropagatorSettings( propagators::body_mass_state ) )
                        prevFinalMasses = missionSegments_[i-1]->getFinalBodyMasses( prevFinalTime );

                // reset the initial time and states of the dynamics simulator
                missionSegments_[i]->setPreliminaryConditions( prevFinalTime, prevFinalStates,
                                                               prevFinalMasses );

            }

            if( nOfDecisionVariables_[i] > 0 )
            {
                newValues.resize( nOfDecisionVariables_[i] );

                for( int j = 0; j < nOfDecisionVariables_[i]; j++)
                {
                    newValues[j] = decisionVariables[accumulator];
                    accumulator++;
                }

                missionSegments_[i]->setPropagationConditions( newValues );
            }

            missionSegments_[i]->dynamicsSimulator_->integrateEquationsOfMotion(
                        missionSegments_[i]->getInitialStates() );

        }

        std::vector< double > returnValue;

        returnValue.push_back( missionSegments_.back()->getObjectiveFunctionValue() );

        return returnValue;
    }

    // Retrieve the boundaries
    //! Retrieve the box-boundaries in PaGMO fashion
    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const
    {
        return boundaries_;
    }

private:

    std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegments_;
    std::pair< std::vector< double >, std::vector< double > > boundaries_;
    std::vector< int > nOfDecisionVariables_;

};



void MissionLinker::optimize( void )
{

    unsigned int i = 0;
    int start_series = 0;


    std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentsBlock;

    do{

        if( start_series == 0 )
        {
            if( i > 0 )
            {
                // get the interpolated end of simulation's epoch
                double prevFinalTime = missionSegmentSettings_[i-1]->getFinalSimulationEpoch();

                // get the interpolated final state of the simulation's epoch
                std::map< std::string, Eigen::Vector6d > prevFinalStates = missionSegmentSettings_[i-1]->
                        getFinalOrbitalStates( prevFinalTime );

                std::map< std::string, double > prevFinalMasses = std::map< std::string, double >();

                if( missionSegmentSettings_[i-1]->hasPropagatorSettings( propagators::body_mass_state ) )
                    prevFinalMasses = missionSegmentSettings_[i-1]->getFinalBodyMasses( prevFinalTime );

                // reset the initial time and states of the dynamics simulator
                missionSegmentSettings_[i]->setPreliminaryConditions( prevFinalTime, prevFinalStates, prevFinalMasses  );
            }

            if( missionSegmentSettings_[i]->decisionVariableSettings_ == NULL )
            {
                missionSegmentSettings_[i]->dynamicsSimulator_->integrateEquationsOfMotion(
                            missionSegmentSettings_[i]->getInitialStates() );

            }
            else{

                missionSegmentsBlock = std::vector< boost::shared_ptr< MissionSegmentSettings > >();

                missionSegmentsBlock.push_back(  missionSegmentSettings_[i] );
                start_series = 1;

                if( missionSegmentSettings_[i]->objectiveFunctionSettings_ != NULL )
                    start_series = 2;
            }
        }
        else if( start_series == 1 )
        {
            missionSegmentsBlock.push_back( missionSegmentSettings_[i] );
            if( missionSegmentSettings_[i]->objectiveFunctionSettings_ != NULL )
                start_series = 2;
        }

        if( start_series == 2 )
        {

            // PAGMO TAKEOVER

            // Create LinkingProblem
            pagmo::problem optimizationProblem = pagmo::problem{
                    LinkingProblem( missionSegmentsBlock ) };
            pagmo::algorithm algo;

            // Set the right algorithm according to the opimization settings
            if( optimizationSettings_->optimizationType_ == global_optimization )
                algo = pagmo::algorithm{ pagmo::de1220() };
            else
                algo = pagmo::algorithm{ pagmo::compass_search() };

            // Set population according to the optimization settings
            pagmo::population::size_type populationSize = optimizationSettings_->populationSize_;
            pagmo::island isl = pagmo::island{ algo, optimizationProblem,
                    populationSize };

            unsigned int counter = 0;

            //Start evolving
            while( true )
            {
                counter ++;
                isl.evolve();
                for( ; isl.status()!=pagmo::evolve_status::idle; )
                    isl.wait();

                if( optimizationSettings_->verbosity_ > 0 )
                {
                    // Show results according to verbosity
                    if( counter % optimizationSettings_->verbosity_ == 0)
                    {
                        std::cout << "Evolution n: " << counter << " Objective Function: " << isl.get_population().champion_f()[0] << "\n";
                        fflush(stdout);
                    }
                }
                // Stop if tolerance has been reached
                if( missionSegmentsBlock.back()->objectiveFunctionSettings_->objectiveValueIsSet_ )
                {
                    if( isl.get_population().champion_f()[0] <= missionSegmentsBlock.back()
                            ->objectiveFunctionSettings_->tolerance_ )
                        break;
                }
                // Stop if max number of evolutions has been reached
                if( optimizationSettings_->stopAtMaxNumberOfEvolutions_ && ( counter == missionSegmentsBlock.back()->objectiveFunctionSettings_->
                                 maxNumberOfEvolutions_ ) )
                        break;
            }

            // END OF PAGMO TAKEOVER


            // Perform the propagation with the optimum values
            unsigned int accumulator = 0;
            std::vector< double > optimizedVariables = isl.get_population().champion_x();
            std::vector< int > nOfDecisionVariables = getNumberOfDecisionVariables( missionSegmentsBlock );

            for( unsigned int j = 0; j < missionSegmentsBlock.size(); j++ )
            {
                if( nOfDecisionVariables[j] > 0 )
                {

                    Eigen::VectorXd newValues;
                    newValues.resize( nOfDecisionVariables[j] );

                    // Partition the optimum values for each missionSegment
                    for( int k = 0; k < nOfDecisionVariables[j]; k++ )
                    {
                        newValues[k] = optimizedVariables[accumulator];
                        accumulator++;
                    }

                    // Set the optimum values as propagation conditions
                    missionSegmentsBlock[j]->decisionVariableValues_ = newValues;
                    missionSegmentsBlock[j]->setPropagationConditions(newValues);

                    // Integrate the mission segment
                    missionSegmentsBlock[j]->dynamicsSimulator_->integrateEquationsOfMotion(
                            missionSegmentsBlock[j]->getInitialStates() );

                }

                // Perform the linkage with the next mission segment
                if( j < missionSegmentsBlock.size() - 1)
                {
                    // Get the final epoch
                    double prevFinalTime = missionSegmentsBlock[j]->getFinalSimulationEpoch();

                    // Get the final states
                    std::map< std::string, Eigen::Vector6d > prevFinalStates = missionSegmentsBlock[j]->
                            getFinalOrbitalStates( prevFinalTime );

                    std::map< std::string, double > prevFinalMasses = std::map< std::string, double >();

                    // Get the final masses
                    if( missionSegmentsBlock[j]->hasPropagatorSettings( propagators::body_mass_state ) )
                        prevFinalMasses = missionSegmentsBlock[j]->getFinalBodyMasses( prevFinalTime );

                    // Reset the initial time and states of the dynamics simulator
                    missionSegmentsBlock[j]->setPreliminaryConditions( prevFinalTime, prevFinalStates,
                                                                       prevFinalMasses  );
                }
            }

            // End of block
            start_series = 0;
        }

        i++;
    }
    while( i < missionSegmentSettings_.size() );
}


}

}
