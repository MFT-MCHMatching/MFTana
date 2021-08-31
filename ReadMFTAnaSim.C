#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "ReadMFTAnaSim.h"
#include "SimulationDataFormat/MCEventHeader.h"

#endif

using MFTAnaSimSATrack = o2::mftana::MFTAnaSimSATrack;
using MFTAnaSimTrack = o2::mftana::MFTAnaSimTrack;
using MFTAnaSimCluster = o2::mftana::MFTAnaSimCluster;
using MFTAnaSimHit = o2::mftana::MFTAnaSimHit;

//_____________________________________________________________________________
void ReadMFTAnaSim(std::string run = "")
{
  gSystem->Load("${ALIBUILD_WORK_DIR}/MFTAna/libMFTAnaSim");

  int verbose = 0;

  double xh, yh, zh, xhi, yhi, zhi, xho, yho, zho;
  double dx, dy, pullx, pully;
  
  // file with histograms
  std::string outfilename = "ReadMFTAnaSim" + run + ".root";
  TFile outFile(outfilename.c_str(),"RECREATE");
  createHistograms();

  // input file (from MFTAna.C)
  std::string inpfilename = "MFTAnaSimTracks" + run + ".root";
  TFile inFile(inpfilename.c_str());
  
  TTree *tree1 = (TTree*)inFile.Get("MFTAnaSimTrack");
  tree1->SetBranchAddress("MFTAnaSimTrack", &anaSimTracksP);
  
  TTree *tree2 = (TTree*)inFile.Get("MFTAnaSimHit");
  tree2->SetBranchAddress("MFTAnaSimHit", &anaSimHitsP);
  
  TTree *tree3 = (TTree*)inFile.Get("MFTAnaSimCluster");
  tree3->SetBranchAddress("MFTAnaSimCluster", &anaSimClustersP);
  
  TTree *tree4 = (TTree*)inFile.Get("MFTAnaSimSATrack");
  tree4->SetBranchAddress("MFTAnaSimSATrack", &anaSimSATracksP);

  tree1->GetEntry(0);
  tree2->GetEntry(0);
  tree3->GetEntry(0);
  tree4->GetEntry(0);

  // MC tracks
  int iMCTrack = 0;
  for (auto& mcTrack : anaSimTracks) {
    if (verbose > 0) {
      printMCTrack(&mcTrack, iMCTrack);
    }
    
    TH1Histos[MCNrOfHits]->Fill(mcTrack.getNHits());
    TH1Histos[MCall_Pt]->Fill(mcTrack.getPt());
    if (mcTrack.getNSATracks() > 0) {
      TH1Histos[MCinSA_Pt]->Fill(mcTrack.getPt());
    }
    if (mcTrack.isTrackable()) {
      TH1Histos[MCtrkbl_Pt]->Fill(mcTrack.getPt());
    }
    if (mcTrack.isFullTrackSA()) {
      auto saTrackIndex = mcTrack.getIntSATrackIndex(0);
      auto& saTrack = anaSimSATracks.at(saTrackIndex);
      const double* iwParam = saTrack.getIwParam();
      const double* owParam = saTrack.getOwParam();
      const double* owXYcov = saTrack.getOwXYcov();
      printf("SAtrack %5d   Iw:   %7.3f  %7.3f  %7.3f   Ow:   %7.3f  %7.3f  %7.3f \n", saTrackIndex, iwParam[0], iwParam[1], iwParam[2], owParam[0], owParam[1], owParam[2]);
      TH1Histos[MCfullSA_Pt]->Fill(mcTrack.getPt());
      int i_cls = 0;
      int ih1 = mcTrack.getFirstHitIndex();
      int ih2 = mcTrack.getLastHitIndex();
      for (int ih = ih1; ih <= ih2; ih++) {
	auto& hit = anaSimHits.at(ih);
	xh = hit.getX();
	yh = hit.getY();
	zh = hit.getZ();
	printf("MCtrack %5d   hit %7d   xyz=  %7.3f  %7.3f  %7.3f \n", iMCTrack, ih, xh, yh, zh);
	if (ih == ih1) {
	  xhi = xh;
	  yhi = yh;
	  zhi = zh;
	}
	if ((ih == ih2) && (mcTrack.getNClusters() == (ih2 - ih1 + 1))) {
	  pullx = (xh - owParam[0]) / owXYcov[0];
	  pully = (yh - owParam[1]) / owXYcov[1];
	  /*
	  if (abs(pullx) > 20. || abs(pully) > 20.) {
	    printf("pullx = %7.4f,   pully = %7.4f !!!\n", pullx, pully);
	  }
	  */
	  TH1Histos[SApullX]->Fill(pullx);
	  TH1Histos[SApullY]->Fill(pully);
	  xho = xh;
	  yho = yh;
	  zho = zh;
	}
	if (mcTrack.getNClusters() == (ih2 - ih1 + 1)) {
	  auto& cluster = anaSimClusters.at(mcTrack.getIntClusIndex(i_cls));
	  dx = cluster.getX() - xh;
	  dy = cluster.getY() - yh;
	  TH2Histos[MCXYdif]->Fill(dx, dy);
	  /*
	  if (abs(dx) > 0.01 || abs(dy) > 0.01) {
	    printf("dx = %7.4f,   dy = %7.3f !!!\n", dx, dy);
	  }
	  */
	  i_cls++;
	}
      } // end hit loop
      if (mcTrack.getNLayers() == 10) {
	TH2Histos[SAXYdifIw]->Fill(xhi - iwParam[0], yhi - iwParam[1]);
	TH2Histos[SAXYdifOw]->Fill(xho - owParam[0], yho - owParam[1]);
      }
    }

    iMCTrack++;
  }

  // SA tracks
  int iSATrack = 0;
  for (auto& saTrack : anaSimSATracks) {
    if (verbose > 0) {
      printSATrack(&saTrack, iSATrack++);
    }
    // calculate pulls
  }
  
  inFile.Close();

  // save histograms
  outFile.cd();
  for (auto& h : TH1Histos) {
    h->Write();
  }
  for (auto& h : TH2Histos) {
    h->Write();
  }
  outFile.Close();
  
}
