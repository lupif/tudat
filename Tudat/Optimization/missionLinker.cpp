#include"missionLinker.h"

#include"pagmo/algorithms/de1220.hpp"
#include"pagmo/island.hpp"
#include"pagmo/problem.hpp"


namespace tudat{

namespace optimization{

// Pagmo problem used ot link the mission segments.
// The purpose of this problem is to retrieve the objective function value at the
// end of a chain of mission segments. The chain is defined as a series of mission
// segments with only one objective function.
struct LinkingProblem{

    LinkingProblem( ){ }

    // Constructor defines the decision variables  boundary vectors and the segmentation associated
    // to each mission segment
    LinkingProblem( std::vector< boost::shared_ptr< MissionSegmentSettings > >& missionSegments ) :
        missionSegments_( missionSegments )
    {
        for( std::vector< boost::shared_ptr< MissionSegmentSettings > >::iterator it = missionSegments_.begin();
             it!= missionSegments_.end(); ++it )
        {
            nOfDecisionVariables_.push_back(0);
            if( *it->decisionVariableSettings_ != NULL )
            {
                std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > singleDecisionVarSettingsVec =
                        *it->decisionVariableSettings_->decisionVariableSettings_;
                for( int i = 0; i<singleDecisionVarSettingsVec.size(); i++ )
                {
                    for( int j = 0; j < singleDecisionVarSettingsVec[i]->boundaries_->
                            getLowerBoundary().rows(); j++ ){
                        boundaries_.first.push_back(
                                    singleDecisionVarSettingsVec[i]->boundaries_->getLowerBoundary()[j] );
                        boundaries_.second.push_back(
                                    singleDecisionVarSettingsVec[i]->boundaries_->getUpperBoundary()[j] );
                        nOfDecisionVariables_.back() += 1;
                    }
                }
            }
        }
    }

    // Retrieve the boundaries
    std::pair< std::vector< double >, std::vector< double > > getBoundaries( void )
    {
        return boundaries_;
    }

    // Calculate the fitness, i.e. the objective function of the last mission segment
    // in the chain
    std::vector< double > fitness( std::vector< double >& decisionVariables )
    {

        double fitnessValue;
        int accumulator = 0;
        Eigen::VectorXd newValues;

        for( int i = 0; i < missionSegments_.size(); i++ ){

            // Set the (interpolated) final states of the previous mission segment
            if( i > 0 ){

                // get the interpolated end of simulation's epoch
                double prevFinalTime = missionSegments_[i-1]->getFinalSimulationEpoch();

                // get the interpolated final state of the simulation's epoch
                std::map< string, Eigen::VectorXd > prevFinalStates = getFinalStates(
                            missionSegments_[i-1]->dynamicsSimulator_, prevFinalTime);

                /*************************************************************************
                 *           INSERT HERE A TRANSLATION BETWEEN REFERENCE FRAMES          *
                 * std::map< string, string > getCentralBodies;                          *
                 * std::map< string, string > getFrameOrientation;                       *
                 * ***********************************************************************/

                // reset the initial time and states of the dynamics simulator
                setPreliminaryConditions( missionSegments_[i]->dynamicsSimulator_,
                                          prevFinalStates, prevFinalTime, );

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
            missionSegments_[i]->dynamicsSimulator_->integrateEquationsOfMotion();

        }

        return missionSegments_.back()->retrieveObjectiveFunctionValue();
    }

private:

    std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegments_;
    std::pair< std::vector< double >, std::vector< double > > boundaries_;
    std::vector< int > nOfDecisionVariables_;

};





}

}
