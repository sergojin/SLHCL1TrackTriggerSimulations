#!/bin/bash

cp ../python/tower/SingleMuonFlatOneOverPt0p0005To0p5_oc_tt25_cfi.py pset.py

for i in $(seq 0 47)
do
# start loop over trigger towers

sed "s@_tower = 25@_tower = ${i}@g" pset.py > ../python/tower/SingleMuonFlatOneOverPt0p0005To0p5_oc_tt${i}_cfi.py

# Call cmsDriver
cmsDriver.py SLHCL1TrackTriggerSimulations/Configuration/python/tower/SingleMuonFlatOneOverPt0p0005To0p5_oc_tt${i}_cfi.py \
    -s GEN,SIM,DIGI:pdigi_valid,L1TrackTrigger \
    --conditions auto:upgradePLS3 \
    --eventcontent RAWSIM \
    --datatier GEN-SIM-DIGI-RAW \
    --beamspot HLLHC \
    --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023TTI,SLHCL1TrackTriggerSimulations/Configuration/customise_pgun.cust_useTrackerOnly \
    --customise_commands "process.VtxSmeared.ZsizeInm = cms.double(0.075 * 1e-6)" \
    --geometry Extended2023TTI \
    --magField 38T_PostLS1 \
    --mc --no_exec --processName RAWSIM \
    -n 1000 --python_filename pset_SingleMuon_oc_tt${i}.py \
    --fileout SingleMuon_oc_tt${i}.root

# Append to the output from cmsDriver
cat <<EOF >> pset_SingleMuon_oc_tt${i}.py



# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Make CRAB recognize the correct output filename
process.RAWSIMoutput.fileName = process.RAWSIMoutput.fileName._value.replace(".root", "_ntuple.root")

# Dump the full python config
with open("dump.py", "w") as f:
    f.write(process.dumpPython())
EOF

# Do not keep the files
if [ ${i} -ne 25 ] && [ ${i} -ne 33 ] && [ ${i} -ne 41 ] ; then
  rm ../python/tower/SingleMuonFlatOneOverPt0p0005To0p5_oc_tt${i}_cfi.py
fi

# end loop over trigger towers
done
