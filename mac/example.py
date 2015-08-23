import sys
from ROOT import gSystem
gSystem.Load("libRecoTool_CCProj")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing CCProj..."

sys.exit(0)

