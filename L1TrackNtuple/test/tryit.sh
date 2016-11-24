#!/bin/bash

ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/ElectronPt2to8_PU0_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > ElectronPt2to8_PU0.txt
#ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/ElectronPt2to8_PU200_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > ElectronPt2to8_PU200.txt
ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/MuonPt2to8_PU0_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > MuonPt2to8_PU0.txt
ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/MuonPt2to8_PU200_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > MuonPt2to8_PU200.txt
ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/PionPt2to8_PU0_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > PionPt2to8_PU0.txt
ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/PionPt2to8_PU200_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > PionPt2to8_PU200.txt
ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/TTbar_PU0_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > TTbar_PU0.txt
ls -v /eos/uscms/store/group/l1upgrades/SLHC/GEN/620_SLHC28p1/blt_projects/NoAnalyzer/TTbar_PU200_part/ntuple_TTI_*.root | sed 's@/eos/uscms/@root://cmseos:1094//@' > TTbar_PU200.txt
