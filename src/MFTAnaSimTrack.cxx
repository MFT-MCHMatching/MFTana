#include "../include/MFTAnaSimTrack.h"

namespace o2::mftana
{

//_____________________________________________________________________________
void MFTAnaSimTrack::init()
{
  for (int i = 0; i < o2::mft::constants::LayersNumber; i++) {
    mLayers[i] = -1;
  }
  for (int i = 0; i < (MCSplitCluster * o2::mft::constants::LayersNumber); i++) {
    mIntSATrackIndex[i] = -1;
    mSATrackMult[i] = 0;
  }
  mNSATracks = 0;
}
  
//_____________________________________________________________________________
bool MFTAnaSimTrack::isFullTrackSA() const
{
  if (mNDisks < 4) {
    return false;
  }
  if (mNSATracks != 1) {
    return false;
  }
  if (mSATrackMult[0] != mNClusters) {
    return false;
  }
  return true;
}
  
}; // end namespace o2::mftana
