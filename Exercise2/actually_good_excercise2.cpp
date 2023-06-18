#include "vector"
#include "omp.h"
#include <immintrin.h>
#include "chrono"
#include "iostream"
#include <fstream>
#include "cmath"
#include <functional>
#include <utility>
#include "math.h"
#include <filesystem>


using namespace std;

vector<double> psi(vector<double> r, double M) {
    vector<double> result(r.size());
    for (int i = 0; i < r.size(); ++i) {
        result[i] = 1 + M/(2*r[i]);
    }
    return result;
}


vector<double> spatial_derivative(double dx, const vector<double> &y, bool neg_bound, bool iszero) {
    int s = y.size();
    int factor;
    int factor2;
    if (neg_bound) {
        factor = -1;
    } else {
        factor = 1;
    }

    if (iszero) {
        factor2 = 0;
    } else {
        factor2 = 1;
    }

    vector<double> result(s);
    for (int i = 2; i < y.size() - 2; ++i) {
        result[i] = (1. / 12 * y[i - 2] - 2./3 * y[i - 1] + 2./3 * y[i + 1] - 1. / 12 * y[i + 2]) / (dx);
    }
    result[0] = (1. / 12 * y[1] * factor - 2./3 * y[0] * factor + 2./3 * y[1] - 1. / 12 * y[2]) / (dx);
    result[1] = (1. / 12 * y[0] * factor - 2./3 * y[0] + 2./3 * y[2] - 1. / 12 * y[3]) / (dx);
    result[s - 1] = (1. / 12 * y[s - 3] - 2./3 * y[s - 2] + 2./3 * factor2 - 1. / 12 * factor2) / (dx);
    result[s - 2] = (1. / 12 * y[s - 4] - 2./3 * y[s - 3] + 2./3 * y[s - 1] - 1. / 12 * factor2) / (dx);

    return result;
}


vector<double> dr_ln_psi(vector<double> r, double M) {
    vector<double> result(r.size());

    for (int i = 0; i < r.size() ; ++i) {
        result[i] = -4*M/(M*r[i]+2*r[i]*r[i]);
    }
    return result;
}

vector<double> ddr_ln_psi(vector<double> r, double M) {
    vector<double> result(r.size());

    for (int i = 0; i < r.size() ; ++i) {
        result[i] = (4*M*(M + 4*r[i]))/(r[i]*r[i]*(M + 2*r[i])*(M + 2*r[i]));
    }
    return result;
}



vector<double> ddt_AB(vector<double> AB, vector<double> alpha, vector<double> K_AB, vector<double> r, double M) {
    vector<double> result(alpha.size());
    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -2 * alpha[i] * AB[i] * K_AB[i];
    }
    return result;
}

vector<double> ddt_D_AB(const vector<double>& D_AB, vector<double> alpha, vector<double> K_AB, vector<double> D_alpha, double dr) {
    vector<double> result(alpha.size());
    vector<double> K_AB_deriv = spatial_derivative(dr, K_AB, false, true);
    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -2 * alpha[i] * (K_AB[i] * D_alpha[i] + K_AB_deriv[i]);
    }
    return result;
}


vector<double> ddt_K_A(vector<double> K_A, vector<double> alpha, vector<double> A, vector<double> K_B,
                       vector<double> D_A, vector<double> D_B, vector<double> D_alpha, vector<double> r, double dr, double M) {
    vector<double> result(alpha.size());
    vector<double> drDalpha(D_alpha.size());
    vector<double> Dalpha_DB_deriv_sum(alpha.size());
    vector<double> drDB(D_alpha.size());
    drDalpha = spatial_derivative(dr, D_alpha, true, true);
    drDB = spatial_derivative(dr, D_B, true, true);

    vector<double> psi_vec = psi(r, M);
    vector<double> dr_ln_psi_vec = dr_ln_psi(r, M);
    vector<double> ddr_ln_psi_vec = ddr_ln_psi(r, M);
    for (int i = 0; i < D_alpha.size(); ++i) {
        Dalpha_DB_deriv_sum[i] = drDalpha[i] + drDB[i];
    }

    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -alpha[i] / (A[i]*pow(psi_vec[i], 4)) * (Dalpha_DB_deriv_sum[i] + ddr_ln_psi_vec[i] + D_alpha[i] * D_alpha[i] - (D_alpha[i] * (D_B[i]+ dr_ln_psi_vec[i])) / 2
                                        + ((D_B[i]+ dr_ln_psi_vec[i]) * (D_B[i]+ dr_ln_psi_vec[i])) / 2 - ((D_A[i]+ dr_ln_psi_vec[i]) * (D_B[i]+ dr_ln_psi_vec[i])) / 2 -
                                        A[i]*pow(psi_vec[i], 4) * K_A[i] * (K_A[i] + 2 * K_B[i]) - 1 / r[i] * ((D_A[i]+ dr_ln_psi_vec[i]) - 2 * (D_B[i]+ dr_ln_psi_vec[i])));
    }
    return result;
}


vector<double> ddt_K_B(vector<double> K_B,vector<double> alpha, vector<double> A, vector<double> B, vector<double> K_A,
                       vector<double> D_A, vector<double> D_B, vector<double> D_alpha, vector<double> r, double dr, double M) {
    vector<double> result(alpha.size());
    vector<double> dr_ln_psi_vec = dr_ln_psi(r, M);
    vector<double> ddr_ln_psi_vec = ddr_ln_psi(r, M);
    vector<double> D_B_deriv = spatial_derivative(dr, D_B, true, true);
    vector<double> psi_vec = psi(r, M);

    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -alpha[i] / (2 * A[i]*pow(psi_vec[i], 4)) * (D_B_deriv[i] + ddr_ln_psi_vec[i] + D_alpha[i] * (D_B[i] + dr_ln_psi_vec[i]) +
                                              (D_B[i] + dr_ln_psi_vec[i]) * (D_B[i] + dr_ln_psi_vec[i]) -
                                              ((D_A[i] + dr_ln_psi_vec[i]) * (D_B[i] + dr_ln_psi_vec[i])) / 2
                                              - 1 / r[i] * ((D_A[i]+ dr_ln_psi_vec[i]) - 2 * D_alpha[i] - 4 * (D_B[i]+ dr_ln_psi_vec[i])) -
                                              (2 * (A[i]*pow(psi_vec[i], 4) - B[i]*pow(psi_vec[i], 4))) / (r[i] * r[i] * B[i]*pow(psi_vec[i], 4)))
                    + alpha[i] * K_B[i] * (K_A[i] + 2 * K_B[i]);
    }
    return result;
}

vector<double> generate_linspace(int N, double dr){
    vector<double> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = (i+1./2)*dr;
    }
    return result;
}


vector<double> MVS(vector<double> v, double k) {
    vector<double> out(v.size());
    for (int i = 0; i < v.size(); i++)
        out[i] = v[i] * k;
    return out;
}


vector<double> addv(vector<double> v, vector<double> w) {
    vector<double> out(v.size());
    for (int i = 0; i < v.size(); i++)
        out[i] = v[i] + w[i];
    return out;
}


vector<vector<double>> RK4(vector<vector<double>> data, double dt, double M, double dr) {
    vector<double> A = data[0];
    vector<double> B = data[1];
    vector<double> D_A = data[2];
    vector<double> D_B = data[3];
    vector<double> K_A = data[4];
    vector<double> K_B = data[5];
    vector<double> alpha = data[6];
    vector<double> D_alpha = data[7];
    vector<double> r = data[8];

    vector<double> w1_A = ddt_AB(A, alpha, K_A, r, M);
    vector<double> w1_B = ddt_AB(B, alpha, K_B, r, M);
    vector<double> w1_D_A = ddt_D_AB(D_A, alpha, K_A, D_alpha, dr);
    vector<double> w1_D_B = ddt_D_AB(D_B, alpha, K_B, D_alpha, dr);
    vector<double> w1_K_A = ddt_K_A(K_A, alpha, A, K_B, D_A, D_B, D_alpha, r, dr, M);
    vector<double> w1_K_B = ddt_K_B(K_B, alpha, A, B, K_A, D_A, D_B, D_alpha, r, dr, M);

    vector<double> w2_A = ddt_AB(addv(A, MVS(w1_A, 0.5 * dt)), alpha, addv(K_A, MVS(w1_K_A, 0.5 * dt)), r, M);
    vector<double> w2_B = ddt_AB(addv(B, MVS(w1_A, 0.5 * dt)), alpha, addv(K_B, MVS(w1_K_B, 0.5 * dt)), r, M);
    vector<double> w2_D_A = ddt_D_AB(addv(D_A, MVS(w1_D_A, 0.5 * dt)), alpha, addv(K_A, MVS(w1_K_A, 0.5 * dt)), D_alpha, dr);
    vector<double> w2_D_B = ddt_D_AB(addv(D_B, MVS(w1_D_B, 0.5 * dt)), alpha, addv(K_B, MVS(w1_K_B, 0.5 * dt)), D_alpha, dr);
    vector<double> w2_K_A = ddt_K_A(addv(K_A, MVS(w1_K_A, 0.5 * dt)), alpha, addv(A, MVS(w1_A, 0.5 * dt)), addv(K_B, MVS(w1_K_B, 0.5 * dt)), addv(D_A, MVS(w1_D_A, 0.5 * dt)), addv(D_B, MVS(w1_D_B, 0.5 * dt)), D_alpha, r, dr, M);
    vector<double> w2_K_B = ddt_K_B(addv(K_B, MVS(w1_K_B, 0.5 * dt)), alpha, addv(A, MVS(w1_A, 0.5 * dt)), addv(B, MVS(w1_B, 0.5 * dt)), addv(K_A, MVS(w1_K_A, 0.5 * dt)), addv(D_A, MVS(w1_D_A, 0.5 * dt)),  addv(D_B, MVS(w1_D_B, 0.5 * dt)), D_alpha, r, dr, M);

    vector<double> w3_A = ddt_AB(addv(A, MVS(w2_A, 0.5 * dt)), alpha, addv(K_A, MVS(w2_K_A, 0.5 * dt)), r, M);
    vector<double> w3_B = ddt_AB(addv(B, MVS(w2_A, 0.5 * dt)), alpha, addv(K_B, MVS(w2_K_B, 0.5 * dt)), r, M);
    vector<double> w3_D_A = ddt_D_AB(addv(D_A, MVS(w2_D_A, 0.5 * dt)), alpha, addv(K_A, MVS(w2_K_A, 0.5 * dt)), D_alpha, dr);
    vector<double> w3_D_B = ddt_D_AB(addv(D_B, MVS(w2_D_B, 0.5 * dt)), alpha, addv(K_B, MVS(w2_K_B, 0.5 * dt)), D_alpha, dr);
    vector<double> w3_K_A = ddt_K_A(addv(K_A, MVS(w2_K_A, 0.5 * dt)), alpha, addv(A, MVS(w2_A, 0.5 * dt)), addv(K_B, MVS(w2_K_B, 0.5 * dt)), addv(D_A, MVS(w2_D_A, 0.5 * dt)), addv(D_B, MVS(w2_D_B, 0.5 * dt)), D_alpha, r, dr, M);
    vector<double> w3_K_B = ddt_K_B(addv(K_B, MVS(w2_K_B, 0.5 * dt)), alpha, addv(A, MVS(w2_A, 0.5 * dt)), addv(B, MVS(w2_B, 0.5 * dt)), addv(K_A, MVS(w2_K_A, 0.5 * dt)), addv(D_A, MVS(w2_D_A, 0.5 * dt)),  addv(D_B, MVS(w2_D_B, 0.5 * dt)), D_alpha, r, dr, M);

    vector<double> w4_A = ddt_AB(addv(A, MVS(w3_A, dt)), alpha, addv(K_A, MVS(w3_K_A, dt)), r, M);
    vector<double> w4_B = ddt_AB(addv(B, MVS(w3_A, dt)), alpha, addv(K_B, MVS(w3_K_B, dt)), r, M);
    vector<double> w4_D_A = ddt_D_AB(addv(D_A, MVS(w3_D_A, dt)), alpha, addv(K_A, MVS(w3_K_A, dt)), D_alpha, dr);
    vector<double> w4_D_B = ddt_D_AB(addv(D_B, MVS(w3_D_B, dt)), alpha, addv(K_B, MVS(w3_K_B, dt)), D_alpha, dr);
    vector<double> w4_K_A = ddt_K_A(addv(K_A, MVS(w3_K_A, dt)), alpha, addv(A, MVS(w3_A, dt)), addv(K_B, MVS(w3_K_B, dt)), addv(D_A, MVS(w3_D_A, dt)), addv(D_B, MVS(w3_D_B, dt)), D_alpha, r, dr, M);
    vector<double> w4_K_B = ddt_K_B(addv(K_B, MVS(w3_K_B, dt)), alpha, addv(A, MVS(w3_A, dt)), addv(B, MVS(w3_B, dt)), addv(K_A, MVS(w3_K_A, dt)), addv(D_A, MVS(w3_D_A, dt)),  addv(D_B, MVS(w3_D_B, dt)), D_alpha, r, dr, M);

    vector<vector<double>> ndata(9, vector<double>(r.size()));

    ndata[0] = MVS(addv(addv(addv(w1_A, MVS(w2_A, 2)), MVS(w3_A, 2)), w4_A), dt / 6);
    ndata[1] = MVS(addv(addv(addv(w1_B, MVS(w2_B, 2)), MVS(w3_B, 2)), w4_B), dt / 6);
    ndata[2] = MVS(addv(addv(addv(w1_D_A, MVS(w2_D_A, 2)), MVS(w3_D_A, 2)), w4_D_A), dt / 6);
    ndata[3] = MVS(addv(addv(addv(w1_D_B, MVS(w2_D_B, 2)), MVS(w3_D_B, 2)), w4_D_B), dt / 6);
    ndata[4] = MVS(addv(addv(addv(w1_K_A, MVS(w2_K_A, 2)), MVS(w3_K_A, 2)), w4_K_A), dt / 6);
    ndata[5] = MVS(addv(addv(addv(w1_K_B, MVS(w2_K_B, 2)), MVS(w3_K_B, 2)), w4_K_B), dt / 6);
    ndata[6] = alpha;
    ndata[7] = D_alpha;
    ndata[8] = r;

    return ndata;
}


// For testing (RK4 stinks)
vector<vector<double>> Euler(vector<vector<double>> data, double dt, double M, double dr) {
    vector<double> A = data[0];
    vector<double> B = data[1];
    vector<double> D_A = data[2];
    vector<double> D_B = data[3];
    vector<double> K_A = data[4];
    vector<double> K_B = data[5];
    vector<double> alpha = data[6];
    vector<double> D_alpha = data[7];
    vector<double> r = data[8];

    vector<double> w1_A = ddt_AB(A, alpha, K_A, r, M);
    vector<double> w1_B = ddt_AB(B, alpha, K_B, r, M);
    vector<double> w1_D_A = ddt_D_AB(D_A, alpha, K_A, D_alpha, dr);
    vector<double> w1_D_B = ddt_D_AB(D_B, alpha, K_B, D_alpha, dr);
    vector<double> w1_K_A = ddt_K_A(K_A, alpha, A, K_B, D_A, D_B, D_alpha, r, dr, M);
    vector<double> w1_K_B = ddt_K_B(K_B, alpha, A, B, K_A, D_A, D_B, D_alpha, r, dr, M);

    vector<vector<double>> ndata(9, vector<double>(r.size()));

    ndata[0] = addv(MVS(w1_A, dt), A);
    ndata[1] = addv(MVS(w1_B, dt), B);
    ndata[2] = addv(MVS(w1_D_A, dt), D_A);
    ndata[3] = addv(MVS(w1_D_B, dt), D_B);
    ndata[4] = addv(MVS(w1_K_A, dt), K_A);
    ndata[5] = addv(MVS(w1_K_B, dt), K_B);
    ndata[6] = alpha;
    ndata[7] = D_alpha;
    ndata[8] = r;

    return ndata;
}



void savetofile(string fname, vector<double> timestepline){
    ofstream of;

    // opening file using ofstream
    of.open(fname, ios::app);
    if (!of){
        ofstream of(fname);
        for (int i = 0; i < timestepline.size()-1; ++i) {
            of << timestepline[i] << ", ";
        }
        of << timestepline[timestepline.size()-1] << "\n";
        of.close();
    }
    else {
        for (int i = 0; i < timestepline.size()-1; ++i) {
            of << timestepline[i] << ", ";
        }
        of << timestepline[timestepline.size()-1] << "\n";
        of.close();
    }
}



void timeevolution(const vector<double>& r, double M, double dt, double N_T, double dr){

    vector<vector<double>> tardata(9, vector<double>(r.size()));

    vector<double> alpha(r.size(), 1);
    vector<double> D_alpha(r.size(), 0);
    vector<double> K_A(r.size(), 0);
    vector<double> K_B(r.size(), 0);
    vector<double> A(r.size(), 1);
    vector<double> B(r.size(), 1);
    vector<double> D_A(r.size(), 0);
    vector<double> D_B(r.size(), 0);

    tardata[0] = A;
    tardata[1] = B;
    tardata[2] = D_A;
    tardata[3] = D_B;
    tardata[4] = K_A;
    tardata[5] = K_B;
    tardata[6] = alpha ;
    tardata[7] = D_alpha;
    tardata[8] = r;
    savetofile("A.csv", tardata[0]);
    savetofile("B.csv", tardata[1]);
    savetofile("D_A.csv", tardata[2]);
    savetofile("D_B.csv", tardata[3]);
    savetofile("K_A.csv", tardata[4]);
    savetofile("K_B.csv", tardata[5]);
    savetofile("alpha.csv", tardata[6]);
    savetofile("D_alpha.csv", tardata[7]);
    savetofile("r.csv", tardata[8]);


    for (int i = 0; i < N_T; ++i) {
        tardata = Euler(tardata, dt, M, dr);
//        tardata = RK4(tardata, dt, M, dr);
        if (i % 200 == 0) {
            savetofile("A.csv", tardata[0]);
            savetofile("B.csv", tardata[1]);
            savetofile("D_A.csv", tardata[2]);
            savetofile("D_B.csv", tardata[3]);
            savetofile("K_A.csv", tardata[4]);
            savetofile("K_B.csv", tardata[5]);
            savetofile("alpha.csv", tardata[6]);
            savetofile("D_alpha.csv", tardata[7]);
            savetofile("r.csv", tardata[8]);
        }
    }

}


int main() {

    std::filesystem::path cwd = std::filesystem::current_path();

    std::vector<std::string> filenames = {
            "A.csv",
            "B.csv",
            "D_A.csv",
            "D_B.csv",
            "K_A.csv",
            "K_B.csv",
            "alpha.csv",
            "D_alpha.csv",
            "r.csv",
    };

    for (const std::string& filename : filenames) {
        std::filesystem::path filepath = cwd / filename;
        if (std::filesystem::exists(filepath)) {
            if (std::filesystem::remove(filepath)) {
                std::cout << "File deleted successfully: " << filepath << std::endl;
            } else {
                std::cerr << "Error deleting the file: " << filepath << std::endl;
                return 1;
            }
        } else {
            std::cout << "File does not exist: " << filepath << std::endl;
        }
    }

    double dr = 0.01;
    vector<double> r = generate_linspace(10000, dr);

    timeevolution(r, 1., 0.005, 2000, dr);

    return 1;
}