#include <vector>
#include <map>

#include "include/MFTAnaSimHit.h"
#include "include/MFTAnaSimCluster.h"
#include "include/MFTAnaSimTrack.h"
#include "include/MFTAnaSimMCTrack.h"
#include "include/MFTAnaSimSATrack.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

enum TH1HistoCodes {
  MCNrOfHits,
  MCall_Pt,
  MCinSA_Pt,
  MCtrkbl_Pt,
  MCfullSA_Pt,
  SApullX,
  SApullY,
  N1Histograms
};

enum TH2HistoCodes {
  SAXYdifIw,
  SAXYdifOw,
  MCXYdif,
  N2Histograms
};

std::map<int, const char *> TH1Names {
  {MCNrOfHits, "MC number of hits"},
  {MCall_Pt, "Pt pf MC tracks with hits"},
  {MCinSA_Pt, "Pt of MC tracks with clusters in SA tracks"},
  {MCtrkbl_Pt, "Pt of MC tracks which are trackable"},
  {MCfullSA_Pt, "Pt of MC tracks which are fully reconstructed"},
  {SApullX, "Pull of x track points"},
  {SApullY, "Pull of y track points"}
};

std::map<int, const char *> TH2Names {
  {SAXYdifIw, "Inward (x,y) diffs to MC hit"},
  {SAXYdifOw, "Outward (x,y) diffs to MC hit"},
  {MCXYdif, "MC (x,y) cls-hit diffs"}    
};

std::map<int, const char *> TH1Titles {
  {MCNrOfHits, "MC number of hits"},
  {MCall_Pt, "Pt pf MC tracks with hits"},
  {MCinSA_Pt, "Pt of MC tracks with clusters in SA tracks"},
  {MCtrkbl_Pt, "Pt of MC tracks which are trackable"},
  {MCfullSA_Pt, "Pt of MC tracks which are fully reconstructed"},
  {SApullX, "Pull of x track points"},
  {SApullY, "Pull of y track points"}
};

std::map<int, const char *> TH2Titles {
  {SAXYdifIw, "Inward (x,y) diffs to MC hit"},
  {SAXYdifOw, "Outward (x,y) diffs to MC hit"},
  {MCXYdif, "MC (x,y) cls-hit diffs"}    
};

std::map<int, std::array<double, 3>> TH1Binning {
  {MCNrOfHits, {50, 0., 50.}},
  {MCall_Pt, {100, 0., 5.}},
  {MCinSA_Pt, {100, 0., 5.}},
  {MCtrkbl_Pt, {100, 0., 5.}},
  {MCfullSA_Pt, {100, 0., 5.}},
  {SApullX, {100, -20., +20.}},
  {SApullY, {100, -20., +20.}}
};

std::map<int, std::array<double, 6>> TH2Binning {
  {SAXYdifIw, {100, -1., +1., 100, -1., +1.}},
  {SAXYdifOw, {100, -1., +1., 100, -1., +1.}},
  {MCXYdif, {100, -0.05, +0.05, 100, -0.05, +0.05}}    
};

std::map<int, const char *> TH1XaxisTitles {
  {MCNrOfHits, "Nr of hits"},
  {MCall_Pt, "pt [GeV/c]"},
  {MCinSA_Pt, "pt [GeV/c]"},
  {MCtrkbl_Pt, "pt [GeV/c]"},
  {MCfullSA_Pt, "pt [GeV/c]"},
  {SApullX, "[cm]"},
  {SApullY, "[cm]"}
};

std::map<int, const char *> TH2XaxisTitles {
  {SAXYdifIw, "dx [cm]"},
  {SAXYdifOw, "dy [cm]"},
  {MCXYdif, "dx [cm]"}
};

std::map<int, const char *> TH2YaxisTitles {
  {SAXYdifIw, "[cm]"},
  {SAXYdifOw, "[cm]"},
  {MCXYdif, "dy [cm]"}
};

std::vector<TH1F*> TH1Histos(N1Histograms);
std::vector<TH2F*> TH2Histos(N2Histograms);

std::vector<o2::mftana::MFTAnaSimTrack> anaSimTracks, *anaSimTracksP = &anaSimTracks;
std::vector<o2::mftana::MFTAnaSimHit> anaSimHits, *anaSimHitsP = &anaSimHits;
std::vector<o2::mftana::MFTAnaSimCluster> anaSimClusters, *anaSimClustersP = &anaSimClusters;
std::vector<o2::mftana::MFTAnaSimSATrack> anaSimSATracks, *anaSimSATracksP = &anaSimSATracks;

//_____________________________________________________________________________
void createHistograms() {  
  auto nHisto = 0;
  for (auto& h : TH1Histos) {
    h = new TH1F(TH1Names[nHisto], TH1Titles[nHisto], (int)TH1Binning[nHisto][0], TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);
    ++nHisto;
  }
  
  auto n2Histo = 0;
  for (auto& h : TH2Histos) {
    h = new TH2F(TH2Names[n2Histo], TH2Titles[n2Histo], (int)TH2Binning[n2Histo][0], TH2Binning[n2Histo][1], TH2Binning[n2Histo][2], (int)TH2Binning[n2Histo][3], TH2Binning[n2Histo][4], TH2Binning[n2Histo][5]);
   h->GetXaxis()->SetTitle(TH2XaxisTitles[n2Histo]);
   h->GetYaxis()->SetTitle(TH2YaxisTitles[n2Histo]);
   //h->GetXaxis()->SetLabelSize(0.05);
   //h->GetXaxis()->SetTitleSize(0.05);
   //h->GetYaxis()->SetLabelSize(0.06);
   //h->GetYaxis()->SetTitleSize(0.06);
   h->SetOption("COLZ");
   ++n2Histo;
  }
}

//_____________________________________________________________________________
void printSATrack(o2::mftana::MFTAnaSimSATrack* saTrack, int i) {
  
  int nPoints, clusEntry;
  auto trkInX = saTrack->getX();
  auto trkInY = saTrack->getY();
  auto trkInZ = saTrack->getZ();
  auto outParam = saTrack->getOutParam();
  auto trkOutX = outParam.getX();
  auto trkOutY = outParam.getY();
  auto trkOutZ = outParam.getZ();
  printf("---------------------------------------------------------------\n");
  printf("SATrack %d points %d disks %d layers %d  x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f  pt  %f\n", i, saTrack->getNPoints(), saTrack->getNDisks(), saTrack->getNLayers(), trkInX, trkInY, trkInZ, trkOutX, trkOutY, trkOutZ, outParam.getPt());
  nPoints = saTrack->getNPoints();
  printf("Number of points %d \n", nPoints);
  if (nPoints > 0) {
    for (int i_p = 0; i_p < nPoints; i_p++) {
      clusEntry = saTrack->getIntClusIndex(i_p);
      auto& cluster = anaSimClusters.at(clusEntry);
      printf("%5d  %5d:   %7.3f   %7.3f   %7.3f   %d \n", i_p, clusEntry, cluster.getX(), cluster.getY(), cluster.getZ(), cluster.getLayer());
    }
  } else {
    printf("... no clusters ...\n");
  }
  
  printf("---------------------------------------------------------------\n");
  for (int iMCTrack = 0; iMCTrack < saTrack->getNMCTracks(); iMCTrack++) {
    printf("MCTrack %d (%d) evn %d \n", saTrack->getIntMCTrackIndex(iMCTrack), saTrack->getMCTrackMult(iMCTrack), saTrack->getMCTrackEvent(iMCTrack));  
  }
  
}

//_____________________________________________________________________________
void printMCTrack(o2::mftana::MFTAnaSimTrack* mcTrack, int i) {

  int nClusters;
  printf("===============================================================\n");
  printf("Event %d Track %d NDisks %d NLayers %d MC %d ", mcTrack->getEvent(), i, mcTrack->getNDisks(), mcTrack->getNLayers(), mcTrack->getMCTrackID());
  printf("layers: ");
  for (int index = 0; index < mcTrack->getNLayers(); index++) {
    printf("%d ", mcTrack->getLayer(index));
  }
  printf("\n");
  
  // hits
  printf("---------------------------------------------------------------\n");
  printf("Number of hits %d , hit range: %d %d \n", mcTrack->getNHits(), mcTrack->getFirstHitIndex(), mcTrack->getLastHitIndex());
  if (mcTrack->getFirstHitIndex() >= 0) {
    for (int iHit = mcTrack->getFirstHitIndex(); iHit <= mcTrack->getLastHitIndex(); iHit++) {
      auto& hit = anaSimHits.at(iHit);
      printf("%5d:   %7.3f   %7.3f   %7.3f \n", iHit, hit.getX(), hit.getY(), hit.getZ());
    }
  }
  
  // clusters
  nClusters = mcTrack->getNClusters();
  printf("---------------------------------------------------------------\n");
  printf("Number of clusters %d \n", nClusters);
  if (nClusters > 0) {
    for (int i_cls = 0; i_cls < nClusters; i_cls++) {
      auto& cluster = anaSimClusters.at(mcTrack->getIntClusIndex(i_cls));
      printf("%5d:   %7.3f   %7.3f   %7.3f \n", i_cls, cluster.getX(), cluster.getY(), cluster.getZ());
    }
  } else {
    printf("... no clusters ...\n");
  }
  
  printf("---------------------------------------------------------------\n");
  for (int iSATrack = 0; iSATrack < mcTrack->getNSATracks(); iSATrack++) {
    printf("SATrack %d (%d) \n", mcTrack->getIntSATrackIndex(iSATrack), mcTrack->getSATrackMult(iSATrack));  
  }

}
