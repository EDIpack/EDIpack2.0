#include <vector>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <edipack2_cbinding.h>  // Include your header file


using namespace std;


// linspace: Linearly spaced values
vector<double> linspace(double start, double end, int num) {
    vector<double> values(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        values[i] = start + i * step;
    }
    return values;
}

// arange: Range with step size
vector<double> arange(double start, double end, double step) {
    vector<double> values;
    for (double val = start; val < end; val += step) {
        values.push_back(val);
    }
    return values;
}

// dens_bethe: Semicircular Density of States
vector<double> dens_bethe(const vector<double>& energies, double W, double de) {
    vector<double> dos(energies.size());
    for (size_t i = 0; i < energies.size(); ++i) {
        double E = energies[i];
        if (abs(E) <= W) {
            dos[i] = de * (2.0 / M_PI / W) * sqrt(1.0 - pow(E / W, 2));
        } else {
            dos[i] = 0.0;
        }
    }
    return dos;
}

// Text colors
string bold_green(const string& s) { return "\033[1;32m" + s + "\033[0m"; }
string bold_yellow(const string& s) { return "\033[1;33m" + s + "\033[0m"; }
string bold_red(const string& s) { return "\033[1;31m" + s + "\033[0m"; }


// Check convergence
bool check_convergence(const complex<double>* Xnew, int dim_Xnew, double eps, int N1, int N2) {
    static complex<double>* Xold = nullptr;
    static int Xold_size = 0;
    static int success = 0;
    static int check = 1;
    bool master = true;

    // Allocate and zero Xold if needed
    if (Xold == nullptr || Xold_size != dim_Xnew) {
        delete[] Xold;
        Xold = new complex<double>[dim_Xnew];
        fill(Xold, Xold + dim_Xnew, complex<double>(0.0, 0.0));
        Xold_size = dim_Xnew;
    }

    double M = 0.0, S = 0.0;
    for (int i = 0; i < dim_Xnew; ++i) {
        M += abs(Xnew[i] - Xold[i]);
        S += abs(Xnew[i]);
    }

    double err = M / S;

    // Update Xold
    memcpy(Xold, Xnew, dim_Xnew * sizeof(complex<double>));

    if (master) {
        ofstream fout("error.err", ios::app);
        fout << setw(5) << check << scientific << setprecision(7) << " " << err << "\n";
    }

    bool convergence = false;
    if (err < eps) {
        ++success;
    } else {
        success = 0;
    }

    if (success > N1) {
        convergence = true;
    }

    if (check >= N2) {
        if (master) {
            ofstream ferr("ERROR.README");
            ferr << "\n";
            cerr << "Not converged after " << N2 << " iterations." << endl;
        }
        convergence = true;
    }

    if (master) {
        if (convergence) {
            cout << bold_green("error=") << scientific << setprecision(7) << err << endl;
        } else if (err < eps) {
            cout << bold_yellow("error=") << scientific << setprecision(7) << err << endl;
        } else {
            cout << bold_red("error=") << scientific << setprecision(7) << err << endl;
        }
    }

    ++check;
    return convergence;
}


int main() {
    char input[] = "inputED.conf"; 
    read_input(input);  
    

    int Le = 1000;
    long int Nb;
    double wmixing = 0.5;
    complex<double> zeta = 0.0;
    double Wband = 1.0;
    double pi = 3.14159265358979323846264338327950288419716939937510;
    bool converged = false;
    int iloop = 0;
    vector<double> Ebands, Dbands;
    double de = 2*Wband/Le;
    complex<double>* wm = new complex<double>[Lmats];
    int64_t d[4] = {Nspin,Nspin,Norb,Norb};
    
    Ebands = linspace(-Wband,Wband,Le);
    Dbands = dens_bethe(Ebands,Wband,de);

    // Initialize frequency arrays
    for (int i = 0; i < Lmats; i++) {
      wm[i] = complex<double>(0.0, pi / ::beta * (2* i + 1));
    }
    
    // Calculate total size of the flattened array
    int total_size = d[0] * d[1] * d[2] * d[3];


    // Allocate a 1D array to hold the flattened Hloc data
    complex<double>* Hloc = new complex<double>[total_size];


    // Call the Fortran subroutine (pass the flattened array and dimensions)
    ed_set_Hloc_single_N4(Hloc, d);

    Nb = get_bath_dimension();
    int64_t bath_dim[1] = {Nb};
    
    double* Bath = new double[Nb];
    double* Bath_prev = new double[Nb];

    int total_size_n5 = total_size * Lmats;
    
    complex<double>* Gimp_mats = new complex<double>[total_size_n5];
    complex<double>* Self_mats = new complex<double>[total_size_n5];
    complex<double>* Gloc_mats = new complex<double>[total_size_n5];
    complex<double>* Delta = new complex<double>[total_size_n5];
    int64_t delta_dim[5] = {Nspin,Nspin,Norb,Norb,Lmats};
 
    for (int i = 0; i < Lmats; ++i) {
          Gimp_mats[i] = complex<double>(0.0,0.0); 
          Self_mats[i] = complex<double>(0.0,0.0); 
          Gloc_mats[i] = complex<double>(0.0,0.0);
          Delta[i] = complex<double>(0.0,0.0);
      }
    
    init_solver_site(Bath, bath_dim);

    while (iloop < Nloop && !converged) {    

      solve_site(Bath, bath_dim, 1, 1);      

      get_gimp_site_n5(Gimp_mats, 0, 0, wm, Lmats, 0);
      get_sigma_site_n5(Self_mats, 0, 0, wm, Lmats, 0);
      

      for (int i = 0; i < Lmats; ++i) {
            zeta = wm[i] + xmu - Self_mats[i];
            Gloc_mats[i] = complex<double>(0.0,0.0);
            for (int j=0; j< Le; j++) {
              Gloc_mats[i] += Dbands[j]/(zeta - Ebands[j]);
            }
        }
        
      for (int i = 0; i < Lmats; ++i) {
        Delta[i] = 0.25 * Wband * Gloc_mats[i];
      }
      
      
      Bath_prev = Bath;
      
      chi2_fitgf_single_normal_n5(Delta, delta_dim, Bath, bath_dim, 1, 0, 1);
      
      for (int i = 0; i <Nb; ++i) {
          Bath[i] = wmixing * Bath[i] + (1.0-wmixing) * Bath_prev[i];
      }
      
      converged = check_convergence(Delta, Lmats, dmft_error, Nsuccess, Nloop);  
    
    }
    
    return 0;
}

