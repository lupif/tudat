/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OPTIMIZATION_MISSION_LINKING_H
#define TUDAT_OPTIMIZATION_MISSION_LINKING_H


#include"missionSegmentSettings.h"


namespace tudat{

namespace optimization{

//! Enum of optimization types
//! \brief The enum of optimization types available in Tudat.
//! (Should be expanded to include multiobjective optimization)
enum OptimizationTypes
{
    global_optimization,
    local_optimization
};


//! Class to define the optimization type, the population size and others
/*!
 * \brief class to define the optimization type, the population size, whether the
 * optimization should end at the maximum number of evolutions defined in the
 * objective function and the verbosity
 */
struct OptimizationSettings
{

    //! Constructor
    /*!
     * \brief Constructor.
     * \param optimizationType If global_optimization the optimization process will use the alghorithm
     * de1220, if local_optimization the compass_search will be used instead
     * \param populationSize Number of individuals in the optimization population
     * \param stopAtMaxNumberOfEvolutions Decide whether the optimization should stop at the maximum number
     * of evolutions defined in the objective function
     * \param verbosity Set to 0 to not receive any message during optimization. n>0 to receive messages every
     * n evolutions
     */
    OptimizationSettings( OptimizationTypes optimizationType = global_optimization,
                          const int populationSize = 32, const bool stopAtMaxNumberOfEvolutions = true,
                          const int verbosity = 0 ) :
        optimizationType_( optimizationType ), populationSize_( populationSize ),
        stopAtMaxNumberOfEvolutions_( stopAtMaxNumberOfEvolutions ), verbosity_( verbosity )
    { }

    //! Storage variable containing the optimization type
    OptimizationTypes optimizationType_;

    //! Storage variable containing the population size
    int populationSize_;

    //! Variable defining whether whether the optimization should stop at the maximum number
    //! of evolutions defined in the objective function
    bool stopAtMaxNumberOfEvolutions_;

    //! Storage variable containing the verbosity
    int verbosity_;

};


//! Class to integrate different mission segmen settings in a chain
struct MissionLinker
{

    //! Constructor
    /*!
     * \brief MissionLinker
     * \param missionSegmentSettings vector of MissionSegmentSettings objects
     * chronologically ordered
     * \param optimizationSettings object to define the optimization type, the population size, whether the
     * optimization should end at the maximum number of evolutions defined in the
     * objective function and the verbosity
     * \param optimizeConfiguration set to true if the user desires to perform
     * the linking and optimization with the costruction
     */
    MissionLinker(  std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentSettings,
                    boost::shared_ptr< OptimizationSettings > optimizationSettings = boost::make_shared< OptimizationSettings >(),
                    const bool optimizeConfiguration = true )
        : missionSegmentSettings_( missionSegmentSettings ), optimizationSettings_(optimizationSettings)
        {

            if( optimizeConfiguration )
                optimize();
        }

    //! Deconstructor
    ~MissionLinker( ){ }

    //! Storage variable containing the vector of mission segment settings
    std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentSettings_;

    //! Object containing information for the optimization process
    boost::shared_ptr< OptimizationSettings > optimizationSettings_;

    //! Function to link and optimize the different mission segments
    /*!
     * Function to link and optimize the different mission segments
     */
    void optimize( void );

private:

};




}

}


#endif
