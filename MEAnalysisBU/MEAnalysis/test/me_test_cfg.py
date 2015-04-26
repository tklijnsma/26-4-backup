import os
from TTH.MEAnalysis.MEAnalysis_cfg import *

process.fwliteInput.pathToFile = cms.string(os.environ["CMSSW_BASE"])

process.fwliteInput.ordering = cms.string("")

process.fwliteInput.samples = cms.VPSet(
    cms.PSet(
        skip     = cms.bool(False),
        name     = cms.string('tthbb_step1'),
        nickName = cms.string('TTH'),
        color    = cms.int32(1),
        xSec     = cms.double(1.0)
    )
)
process.fwliteInput.outFileName = cms.string("test_step2.root")
process.fwliteInput.debug = cms.untracked.int32(3)
process.fwliteInput.printout = cms.untracked.int32(1)
process.fwliteInput.evLimits = cms.vint32(0, 500)
process.fwliteInput.ntuplizeAll = cms.untracked.int32(1)
process.fwliteInput.cutLeptons = cms.untracked.bool(False)
process.fwliteInput.cutJets = cms.untracked.bool(False)
process.fwliteInput.cutWMass = cms.untracked.bool(False)
process.fwliteInput.cutBTagShape = cms.untracked.bool(False)
