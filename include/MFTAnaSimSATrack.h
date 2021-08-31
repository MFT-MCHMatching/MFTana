#ifndef MFT_ANA_SA_TRACK
#define MFT_ANA_SA_TRACK

#include "MFTBase/Constants.h"
#include "DataFormatsMFT/TrackMFT.h"

namespace o2::mftana
{

constexpr int SASplitCluster = 4;   ///< Maximum split of a cluster with the same MC track ID

class MFTAnaSimSATrack : public o2::mft::TrackMFT
{
 public:
  MFTAnaSimSATrack() = default;
  ~MFTAnaSimSATrack() = default;
  MFTAnaSimSATrack& operator=(const MFTAnaSimSATrack&) = default;

  void init();
  
  void setNDisks(int nd) { mNDisks = nd; }
  int getNDisks() const { return mNDisks; }
  
  void setNLayers(int nl) { mNLayers = nl; }
  int getNLayers() const { return mNLayers; }
  
  void setNPoints(int np) { mNPoints = np; }
  int getNPoints() const { return mNPoints; }
  
  void setLayer(int ipoint, int layer) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mLayers[ipoint] = layer;
  }
  int getLayer(int ipoint) const { return mLayers[ipoint]; }

  void setIntClusIndex(int ipoint, int index) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mIntClusIndex[ipoint] = index;
  }
  int getIntClusIndex(int ipoint) const {
    assert(ipoint < o2::mft::constants::LayersNumber);
    return mIntClusIndex[ipoint];
  }
  
  void setEvent(int ipoint, int evnID) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mEventID[ipoint] = evnID;
  }
  int getEvent(int ipoint) const {
    assert(ipoint < o2::mft::constants::LayersNumber);
    return mEventID[ipoint];
  }

  void setMCTrackID(int ipoint, int trkID) {
    assert(ipoint < o2::mft::constants::LayersNumber);
    mMCTrackID[ipoint] = trkID;
  }
  int getMCTrackID(int ipoint) const {
    assert(ipoint < o2::mft::constants::LayersNumber);
    return mMCTrackID[ipoint];
  }

  void copy(o2::mft::TrackMFT& track);
  
  int getNMCTracks() const { return mNMCTracks; }
  
  void addMCTrack(int evn, int mcid);
  
  void setIntMCTrackIndex(int i, int index);
  
  int getMCTrackEvent(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mMCTrackEvent[i];
  }
  
  int getMCTrackIndex(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mMCTrackIndex[i];
  }
  
  int getMCTrackMult(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mMCTrackMult[i];
  }

  int getIntMCTrackIndex(int i) const {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    return mIntMCTrackIndex[i];
  }

  void setIwParam(double parx, double pary, double parz) { mIwParam[0] = parx; mIwParam[1] = pary; mIwParam[2] = parz; }
  double* getIwParam() { return mIwParam; }
  
  void setOwParam(double parx, double pary, double parz) { mOwParam[0] = parx; mOwParam[1] = pary; mOwParam[2] = parz; }
  double* getOwParam() { return mOwParam; }
  
  void setIwXYcov(double covx, double covy) { mIwXYcov[0] = covx; mIwXYcov[1] = covy; }
  double* getIwXYcov() { return mIwXYcov; }
  
  void setOwXYcov(double covx, double covy) { mOwXYcov[0] = covx; mOwXYcov[1] = covy; }
  double* getOwXYcov() { return mOwXYcov; }
  
  void cookTrack();
  
 private: 
  int mNDisks = 0;   ///< Number of MFT disks
  int mNLayers = 0;   ///< Number of MFT layers
  int mNPoints = 0;   ///< Number of points in the track
  int mLayers[o2::mft::constants::LayersNumber];   ///< ID of the layers
  int mIntClusIndex[o2::mft::constants::LayersNumber];   ///< Internal index for the attached clusters
  int mEventID[o2::mft::constants::LayersNumber];   ///< ID of the events to which the clusters belong
  int mMCTrackID[o2::mft::constants::LayersNumber];   ///< ID of the MC tracks which contribute to the points
  
  int mMCTrackEvent[SASplitCluster * o2::mft::constants::LayersNumber];   ///< List of the event IDs which contribute with MC tracks
  int mMCTrackIndex[SASplitCluster * o2::mft::constants::LayersNumber];   ///< List of MC track IDs which contribute with clusters to this SA track
  int mMCTrackMult[SASplitCluster * o2::mft::constants::LayersNumber];   ///< Multiplicity of MC track IDs which contribute with clusters to this SA track
  int mIntMCTrackIndex[SASplitCluster * o2::mft::constants::LayersNumber];   ///< List of MC track (internal) tree indexes which contribute with clusters to this SA track
  int mNMCTracks = 0;   ///< Number of MC tracks which contribute with clusters to this SA track
  double mIwParam[3];   ///< Track parameters "inward" (vertex)
  double mOwParam[3];   ///< Track parameters "outward" (MCH)
  double mIwXYcov[2];   ///< Track x, y covariances "inward" (vertex)
  double mOwXYcov[2];   ///< Track x, y covariances "outward" (MCH)
};

//_____________________________________________________________________________
inline void MFTAnaSimSATrack::setIntMCTrackIndex(int i, int index)
{
  mIntMCTrackIndex[i] = index;
}

//_____________________________________________________________________________
inline void MFTAnaSimSATrack::addMCTrack(int evn, int mcid)
{
  int i;
  for (i = 0; i < mNMCTracks; i++) {
    assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
    if (mMCTrackEvent[i] == evn && mMCTrackIndex[i] == mcid ) {
      mMCTrackMult[i]++;
      return;
    }
  }
  //printf("%d  %d  %d \n", i, evn, mcid);
  assert(i < (SASplitCluster * o2::mft::constants::LayersNumber));
  mMCTrackEvent[i] = evn;
  mMCTrackIndex[i] = mcid;
  mMCTrackMult[i]++;
  mNMCTracks++;
}

};

#endif // MFT_ANA_SA_TRACK
