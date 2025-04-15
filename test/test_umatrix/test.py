import numpy as np
from edipy2 import global_env as ed




#READ ED INPUT:
ed.read_input("inputED.conf")


#Generate hk and hloc
Hloc=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb),dtype=complex)


#SETUP SOLVER
ed.set_hloc(Hloc)
Nb=ed.get_bath_dimension()
bath = ed.init_solver()


