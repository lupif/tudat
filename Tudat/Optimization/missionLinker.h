#ifndef TUDAT_OPTIMIZATION_MISSION_LINKING_H
#define TUDAT_OPTIMIZATION_MISSION_LINKING_H


#include"missionSegmentSettings.h"


namespace tudat{

namespace optimization{

struct MissionLinker
{

    MissionLinker(  std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentSettings  )
        : missionSegmentSettings_( missionSegmentSettings )
        { }

    std::vector< std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentSettings_;

private:








};


void optimize( start, finish );




}

}


#endif
