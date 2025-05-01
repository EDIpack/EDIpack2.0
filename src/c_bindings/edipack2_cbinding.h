#ifndef EDIPACK2_CBINDING
#define EDIPACK2_CBINDING

#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
*This flag is `1` if real-space DMFT support enabled, `0` othervise
*/ 
extern int has_ineq;

/*!
\rst
Inverse temperature
Interfaces to :f:var:`beta`.
\endrst
*/ 
extern double beta;

/*!
\rst
Convergence threshold
Interfaces to :f:var:`dmft_error`.
\endrst
*/ 
extern double dmft_error;

/*!
\rst
Twin sector simmetry flag
Interfaces to :f:var:`ed_twin`.
\endrst
*/ 
extern bool ed_twin;

/*!
\rst
Flag setting orbital-resolved sectorization of Fock space Hamiltoniaj
Interfaces to :f:var:`ed_total_ud`.
\endrst
*/ 
extern bool ed_total_ud;

/*!
\rst
Real axis broadening
Interfaces to :f:var:`eps`.
\endrst
*/ 
extern double eps;

/*!
\rst
Hund's coupling.
Interfaces to :f:var:`Jh`.
\endrst
*/ 
extern double Jh;

/*!
\rst
Pair-hopping interaction strength.
Interfaces to :f:var:`Jp`.
\endrst
*/ 
extern double Jp;

/*!
\rst
Spin-exchange interaction strength.
Interfaces to :f:var:`Jx`.
\endrst
*/ 
extern double Jx;

/*!
\rst
Number of frequencies used in :math:`\chi^{2}` fit
Interfaces to :f:var:`Lfit`.
\endrst
*/ 
extern int Lfit;

/*!
\rst
Number of Matsubara frequencies.
Interfaces to :f:var:`Lmats`.
\endrst
*/ 
extern int Lmats;

/*!
\rst
Number of real frequencies.
Interfaces to :f:var:`Lreal`.
\endrst
*/ 
extern int Lreal;

/*!
\rst
Interfaces to :f:var:`Lpos`.
\endrst
*/ 
extern int Lpos;

/*!
\rst
Number of imaginary-time points.
Interfaces to :f:var:`Ltau`.
\endrst
*/ 
extern int Ltau;

/*!
\rst
Number of bath sites
Interfaces to :f:var:`Nbath`.
\endrst
*/ 
extern int Nbath;

/*!
\rst
Number of DMFT loops
Interfaces to :f:var:`Nloop`.
\endrst
*/ 
extern int Nloop;

/*!
\rst
Number of orbitals
Interfaces to :f:var:`Norb`.
\endrst
*/ 
extern int Norb;

/*!
\rst
Number of phonons.
Interfaces to :f:var:`Nph`.
\endrst
*/ 
extern int Nph;

/*!
\rst
Occupation value for fixed-density calculations
Interfaces to :f:var:`nread`.
\endrst
*/ 
extern double nread;

/*!
\rst
Number of spins
Interfaces to :f:var:`Nspin`.
\endrst
*/ 
extern int Nspin;

/*!
\rst
Number of iterations under convergence threshod
Interfaces to :f:var:`Nsuccess`.
\endrst
*/ 
extern int Nsuccess;

/*!
\rst
Symmetry-breaking field.
Interfaces to :f:var:`sb_field`.
\endrst
*/ 
extern double sb_field;

/*!
\rst
Hubbard interaction strength.
Interfaces to :f:var:`Uloc`.
\endrst
*/ 
extern double Uloc[5];

/*!
\rst
Density-density Kanamori interaction strength.
Interfaces to :f:var:`Ust`.
\endrst
*/ 
extern double Ust;

/*!
\rst
Upper bound on real axis frequencies
Interfaces to :f:var:`wfin`.
\endrst
*/ 
extern double wfin;

/*!
\rst
Upper bound on real axis frequencies
Interfaces to :f:var:`wini`.
\endrst
*/ 
extern double wini;

/*!
\rst
Interfaces to :f:var:`xmax`.
\endrst
*/ 
extern double xmax;

/*!
\rst
Interfaces to :f:var:`xmin`.
\endrst
*/ 
extern double xmin;


/*!
\rst
Chemical potential
Interfaces to :f:var:`xmu`.
\endrst
*/ 
extern double xmu;




/*!
\rst
This function reads the input. Interfaces to :f:func:`ed_read_input`
\endrst
*
* @param instr: the input file name
*/
void read_input(char *instr);


/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/ed_aux_funx/ed_set_Hloc`.
Rank-2 array variant for single-site DMFT
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
*/
void ed_set_Hloc_single_N2(std::complex<double> *Hloc, 
                           int64_t *d
                           );

/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/ed_aux_funx/ed_set_Hloc`.
Rank-4 array variant for single-site DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
*/       
void ed_set_Hloc_single_N4(std::complex<double> *Hloc, 
                           int64_t *d
                           );
                           
                           
/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/e2i_aux_funx/ed_set_Hloc`.
Rank-2 array variant for real-space DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
* @param Nlat: number of inequivalent sites
*/                      
void ed_set_Hloc_lattice_N2(std::complex<double> *Hloc,
                            int64_t *d,
                            int Nlat
                            );

/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/e2i_aux_funx/ed_set_Hloc`.
Rank-2 array variant for real-space DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
* @param Nlat: number of inequivalent sites
*/ 
void ed_set_Hloc_lattice_N3(std::complex<double> *Hloc,
                            int64_t *d,
                            int Nlat
                            );

/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/e2i_aux_funx/ed_set_Hloc`.
Rank-2 array variant for real-space DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
* @param Nlat: number of inequivalent sites
*/ 
void ed_set_Hloc_lattice_N5(std::complex<double> *Hloc, 
                            int64_t *d, 
                            int Nlat
                            );



/*!
\rst
This function gets the dimension of the user-accessible bath array.
Interfaces to :f:func:`get_bath_dimension`.
\endrst
*
* @return the dimension of the bath array
*/ 
int  get_bath_dimension(void);

/*!
\rst
This function updates a variable (usually the chemical potential) trying to 
achieve a desired density set by :f:var:`nread`.
Interfaces to :f:func:`ed_search-variable`.
\endrst
*
* @param var: the variable to be adjusted
* @param ntmp: the density value at a given iteration
* @param converged: the DMFT loop convergence status.
*/ 
void search_variable(double *var,
                    double *ntmp, 
                    int64_t *converged
                    );


/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hreplica`.
Rank-5 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/ 
void init_Hreplica_symmetries_d5(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );
                                 
/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hreplica`.
Rank-3 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                                 
void init_Hreplica_symmetries_d3(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );
                                 
/*!
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hgeneral`.
Rank-5 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                             
void init_Hgeneral_symmetries_d5(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );

/*!
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hgeneral`.
Rank-3 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                               
void init_Hgeneral_symmetries_d3(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );


/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hreplica`.
Rank-5 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/   
void init_Hreplica_symmetries_lattice_d5(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );
                                         
/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hreplica`.
Rank-3 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                                            
void init_Hreplica_symmetries_lattice_d3(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );
                                         
/*!
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hgeneral`.
Rank-5 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/ 
void init_Hgeneral_symmetries_lattice_d5(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );
/*!                                         
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hgeneral`.
Rank-5 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the matrix linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                                       
void init_Hgeneral_symmetries_lattice_d3(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );

void break_symmetry_bath_site(double *bath, 
                              int64_t *dim_bath, 
                              double field, 
                              double sgn, 
                              int sav
                              );
void break_symmetry_bath_ineq(double *bath, 
                              int64_t *dim_bath, 
                              double field, 
                              double sgn, 
                              int sav
                              );

void spin_symmetrize_bath_site(double *bath, 
                               int64_t *dim_bath, 
                               int sav
                               );
void spin_symmetrize_bath_ineq(double *bath, 
                               int64_t *dim_bath, 
                               int sav
                               );

void orb_symmetrize_bath_site(double *bath, 
                              int64_t *dim_bath, 
                              int orb1, 
                              int orb2, 
                              int sav
                              );
void orb_symmetrize_bath_ineq(double *bath, 
                              int64_t *dim_bath, 
                              int orb1, 
                              int orb2, 
                              int sav
                              );

void orb_equality_bath_site(double *bath, 
                            int64_t *dim_bath, 
                            int indx, 
                            int sav
                            );
void orb_equality_bath_ineq(double *bath, 
                            int64_t *dim_bath, 
                            int indx, 
                            int sav
                            );

void ph_symmetrize_bath_site(double *bath, 
                             int64_t *dim_bath, 
                             int sav
                             );
void ph_symmetrize_bath_ineq(double *bath, 
                             int64_t *dim_bath, 
                             int sav
                             );

void save_array_as_bath_site(double *bath, 
                             int64_t *dim_bath
                             );
void save_array_as_bath_ineq(double *bath, 
                             int64_t *dim_bath
                             );

void chi2_fitgf_single_normal_n3(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );
void chi2_fitgf_single_normal_n5(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );
                                
void chi2_fitgf_single_superc_n3(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );
void chi2_fitgf_single_superc_n5(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );

void chi2_fitgf_lattice_normal_n3(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );
void chi2_fitgf_lattice_normal_n4(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );
void chi2_fitgf_lattice_normal_n6(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );

void chi2_fitgf_lattice_superc_n3(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );
void chi2_fitgf_lattice_superc_n4(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );
void chi2_fitgf_lattice_superc_n6(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );
                                 
                                 
void ed_get_dens_n1(std::complex<double> *self);
void ed_get_dens_n2(std::complex<double> *self, int Nlat);

void ed_get_mag_n2(std::complex<double> *self);
void ed_get_mag_n3(std::complex<double> *self, int Nlat);

void ed_get_docc_n1(std::complex<double> *self);
void ed_get_docc_n2(std::complex<double> *self, int Nlat);

void ed_get_phisc_n2(std::complex<double> *self);
void ed_get_phisc_n3(std::complex<double> *self, int Nlat);

void ed_get_eimp_n1(std::complex<double> *self);
void ed_get_eimp_n2(std::complex<double> *self, int Nlat);


void get_sigma_site_n3(std::complex<double> *self, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );
void get_sigma_site_n5(std::complex<double> *self, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );
                       
void get_sigma_lattice_n3(std::complex<double> *self, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
void get_sigma_lattice_n4(std::complex<double> *self, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
void get_sigma_lattice_n6(std::complex<double> *self, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );

void get_gimp_site_n3(std::complex<double> *gimp, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );
void get_gimp_site_n5(std::complex<double> *gimp, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );
                       
void get_gimp_lattice_n3(std::complex<double> *gimp, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
void get_gimp_lattice_n4(std::complex<double> *gimp, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
void get_gimp_lattice_n6(std::complex<double> *gimp, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );

void get_g0and_n3(std::complex<double> *warray, 
                  int64_t *dim_warray, 
                  double *bath, 
                  int *dim_bath, 
                  std::complex<double> *G0and, 
                  int64_t *dim_g0and, 
                  char *axis, 
                  char *typ
                  );
/*!
\rst
This function returns the Weiss field in the single-impurity case.
Interfaces to :f:func:`f/ed_io/ed_get_g0and` for a rank-5 matrix.
\endrst
* @param warray: Array of frequencies
* @param dim_warray: Dimension of the array of frequencies
* @param bath: User-accessible bath array
* @param dim_bath: Dimension of the bath array
* @param G0and: the Weiss field matrix
* @param dim_g0and: the dimensions of the Weiss field matrix
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
*
*/                
void get_g0and_n5(std::complex<double> *warray,
                  int64_t *dim_warray, 
                  double *bath, 
                  int *dim_bath, 
                  std::complex<double> *G0and, 
                  int64_t *dim_g0and, 
                  char *axis = (char *)"m", 
                  char *typ = (char *)"n"
                  );
                  
void get_delta_n3(std::complex<double> *warray, 
                  int64_t *dim_warray, 
                  double *bath, 
                  int *dim_bath, 
                  std::complex<double> *Delta, 
                  int64_t *dim_delta, 
                  char *axis, 
                  char *typ
                  );
void get_delta_n5(std::complex<double> *warray, 
                  int64_t *dim_warray, 
                  double *bath, 
                  int *dim_bath, 
                  std::complex<double> *Delta, 
                  int64_t *dim_delta, 
                  char *axis, 
                  char *typ
                  );

void ed_get_spinchi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );
void ed_get_denschi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );
void ed_get_pairchi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );
void ed_get_exctchi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );

void init_solver_site(double *bath, int64_t *dim_bath);
void init_solver_ineq(double *bath, int64_t *dim_bath);

void solve_site(double *bath, 
                int64_t *dim_bath, 
                int flag_gf, 
                int flag_mpi
                );
void solve_ineq(double *bath, 
                int64_t *dim_bath, 
                int mpi_lanc, 
                int flag_gf
                );

/*! 
\cond
void finalize_solver (@WITH_INEQ_HEADER@)
\endcond
*/

void reset_umatrix();

void add_twobody_operator(int o1,
                          int s1,
                          int o2,
                          int s2,
                          int o3,
                          int s3,
                          int o4,
                          int s4,
                          double U
                     );

#ifdef __cplusplus
}
#endif
#endif
