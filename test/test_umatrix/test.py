import numpy as np
from edipy2 import global_env as ed
import sys



#READ ED INPUT:
if(sys.argv[1]):
    ed.read_input(sys.argv[1])
else:
    ed.read_input("inputED.conf")


#Generate hk and hloc
Hloc=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb),dtype=complex)


#SETUP SOLVER
ed.set_hloc(Hloc)
Nb=ed.get_bath_dimension()
bath = ed.init_solver()

ed.solve(bath)


