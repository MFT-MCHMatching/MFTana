#include "../include/MFTAnaSimSATrack.h"

namespace o2::mftana
{

//_____________________________________________________________________________
void MFTAnaSimSATrack::init()
{
  for (int i = 0; i < o2::mft::constants::LayersNumber; i++) {
    mLayers[i] = -1;
  }
  for (int i = 0; i < (SASplitCluster * o2::mft::constants::LayersNumber); i++) {
    mMCTrackEvent[i] = -1;
    mMCTrackIndex[i] = -1;
    mMCTrackMult[i] = 0;
    mIntMCTrackIndex[i] = -1;
  }
  mNMCTracks = 0;
}
  
//_____________________________________________________________________________
void MFTAnaSimSATrack::copy(o2::mft::TrackMFT& track)
{
  this->setOutParam(track.getOutParam());
  this->setX(track.getX());
  this->setY(track.getY());
  this->setZ(track.getZ());
  this->setCharge(track.getCharge());
  this->setTrackChi2(track.getTrackChi2());
}
  
//_____________________________________________________________________________
void MFTAnaSimSATrack::cookTrack()
{
  //printf("cookTrack\n");
  for (int ip = 0; ip < mNPoints; ip++) {
    addMCTrack(getEvent(ip), getMCTrackID(ip));
  }
}

}; // end namespace o2::mftana
