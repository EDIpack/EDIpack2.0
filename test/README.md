### TEST

We test all cases using 2-orbital problems:

* `BATH_TYPE = NORMAL`
    * `ED_MODE = NORMAL`
    * `ED_MODE = SUPERC` including order parameter $\langle c_{a\uparrow}c_{a\downarrow} \rangle$.
    * `ED_MODE = NONSU2` including in-plane magnetization  $S_x$.

* `BATH_TYPE = HYBRID`
    * `ED_MODE = NORMAL`
    * `ED_MODE = SUPERC` including order parameter $\langle c_{a\uparrow}c_{a\downarrow} \rangle$.
    * `ED_MODE = NONSU2` including in-plane magnetization  $S_x$.

* `BATH_TYPE = REPLICA`
    * `ED_MODE = NORMAL`
    * `ED_MODE = SUPERC` including intra and inter-orbital SC order parameter $\langle c_{a\uparrow}c_{b\downarrow} \rangle$.
    * `ED_MODE = NONSU2` icluding excitonic terms $S_0,T_{x,y,z}$.    

* `BATH_TYPE = GENERAL`
    * `ED_MODE = NORMAL`
    * `ED_MODE = SUPERC` including intra and inter-orbital SC order parameter $\langle c_{a\uparrow}c_{b\downarrow} \rangle$.
    * `ED_MODE = NONSU2` icluding excitonic terms $S_0,T_{x,y,z}$.    


Each code includes tests for different options, according to the following options:  
* `ED_SPARSE_H     = T/F` *(sparse or direct solver)*  
* `ED_READ_UMATRIX = T/F` *(interaction controlled by local parameter* `Uloc`,`Ust`,`Jh`,`Jx`,`Jp`  *or read from file)*
* `ED_USE_KANAMORI = T/F` *(either read from file or use* `add_twobody_operator`   *procedure to set the interaction terms)*

###BONUS  
We also provide a quick test of the `EDIpack2ineq` functionalities for a two impurities case, corresponding to AFM doubled unit cell with two atoms. 

###Remark
The are currently no tests with phonons
