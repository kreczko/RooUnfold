from ROOT import gSystem
gSystem.Load('libRooUnfold.so')
from ROOT import RooUnfold, RooUnfoldSvd, RooUnfoldResponse
from ROOT import RooUnfoldBayes, RooUnfoldBinByBin, RooUnfoldErrors
from ROOT import RooUnfoldInvert, RooUnfoldParms, RooUnfoldTUnfold
from ROOT import TauSVDUnfold 
from ROOT import RooUnfoldUtils
assert  'calculate_tau_scan_points' in dir(RooUnfoldUtils)