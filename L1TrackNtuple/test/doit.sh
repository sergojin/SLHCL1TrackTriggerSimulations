#!/bin/bash

mkdir -p TrkPlots/ TrkPlots2/

# L1TrackNtuplePlot(TString type, int TP_select_injet=0, int TP_select_pdgid=0, int TP_select_eventid=0, float TP_minPt=3.0, float TP_maxPt=100.0, float TP_maxEta=2.4)

root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt2to8_PU0.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt2to8_PU140.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt2to8_PU200.txt\",0,11\)
root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt8to100_PU0.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt8to100_PU140.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt8to100_PU200.txt\",0,11\)
root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt2to8_PU0.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt2to8_PU140.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt2to8_PU200.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt8to100_PU0.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt8to100_PU140.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt8to100_PU200.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt2to8_PU0.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt2to8_PU140.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt2to8_PU200.txt\",0,211\)
root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt8to100_PU0.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt8to100_PU140.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt8to100_PU200.txt\",0,211\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",0,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",0,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",0,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",1,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",1,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",1,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",2,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",2,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",2,0\)

sed -i 's@#define AMTTNROADS 99999999@#define AMTTNROADS 200@' L1TrackNtuplePlot.C

root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt2to8_PU0.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt2to8_PU140.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt2to8_PU200.txt\",0,11\)
root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt8to100_PU0.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt8to100_PU140.txt\",0,11\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"ElectronPt8to100_PU200.txt\",0,11\)
root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt2to8_PU0.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt2to8_PU140.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt2to8_PU200.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt8to100_PU0.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt8to100_PU140.txt\",0,13\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"MuonPt8to100_PU200.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt2to8_PU0.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt2to8_PU140.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt2to8_PU200.txt\",0,211\)
root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt8to100_PU0.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt8to100_PU140.txt\",0,211\)
#root -l -b -q L1TrackNtuplePlot.C+\(\"PionPt8to100_PU200.txt\",0,211\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",0,13\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",0,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",0,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",0,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",1,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",1,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",1,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU0.txt\",2,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU140.txt\",2,0\)
root -l -b -q L1TrackNtuplePlot.C+\(\"TTbar_PU200.txt\",2,0\)

sed -i 's@#define AMTTNROADS 200@#define AMTTNROADS 99999999@' L1TrackNtuplePlot.C

root -l -b -q PlotL1iso.C+\(\"TTbar_PU0.txt\"\)
root -l -b -q PlotL1iso.C+\(\"TTbar_PU140.txt\"\)
root -l -b -q PlotL1iso.C+\(\"TTbar_PU200.txt\"\)
root -l -b -q PlotL1PV.C+\(\"TTbar_PU0.txt\"\)
root -l -b -q PlotL1PV.C+\(\"TTbar_PU140.txt\"\)
root -l -b -q PlotL1PV.C+\(\"TTbar_PU200.txt\"\)
cp output_L1Iso_TTbar_PU0_WithoutTruncation.root output_L1Iso_TTbar_PU0_WithTruncation.root
cp output_L1Iso_TTbar_PU140_WithoutTruncation.root output_L1Iso_TTbar_PU140_WithTruncation.root
cp output_L1Iso_TTbar_PU200_WithoutTruncation.root output_L1Iso_TTbar_PU200_WithTruncation.root
cp output_L1PV_TTbar_PU0_WithoutTruncation.root output_L1PV_TTbar_PU0_WithTruncation.root
cp output_L1PV_TTbar_PU140_WithoutTruncation.root output_L1PV_TTbar_PU140_WithTruncation.root
cp output_L1PV_TTbar_PU200_WithoutTruncation.root output_L1PV_TTbar_PU200_WithTruncation.root

rm -f TrkPlots/*.eps TrkPlots2/*.eps
