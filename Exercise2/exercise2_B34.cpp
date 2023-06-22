#include "vector"
#include "chrono"
#include "iostream"
#include <fstream>
#include "cmath"
#include <filesystem>
#include <iomanip>


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
    s -= 1;
    result[0] = (1. / 12 * y[1] * factor - 2./3 * y[0] * factor + 2./3 * y[1] - 1. / 12 * y[2]) / (dx);
    result[1] = (1. / 12 * y[0] * factor - 2./3 * y[0] + 2./3 * y[2] - 1. / 12 * y[3]) / (dx);
    result[s] = (1. / 12 * y[s - 2] - 2./3 * y[s - 1] + 2./3 * factor2 - 1. / 12 * factor2) / (dx);
    result[s - 1] = (1. / 12 * y[s - 3] - 2./3 * y[s - 2] + 2./3 * y[s] - 1. / 12 * factor2) / (dx);

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


vector<double> vec_log(vector<double> vec){
    vector<double> result(vec.size());
    for (int i = 0; i < vec.size() ; ++i) {
        result[i] = log(vec[i]);
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

vector<double> ddt_D_AB(const vector<double>& D_AB, vector<double> alpha, vector<double> K_AB, double dr) {
    vector<double> result(alpha.size());
    vector<double> D_alpha = spatial_derivative(dr, vec_log(alpha), true, true);
    vector<double> K_AB_deriv = spatial_derivative(dr, K_AB, false, true);
    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -2 * alpha[i] * (K_AB[i] * D_alpha[i] + K_AB_deriv[i]);
    }
    return result;
}


vector<double> ddt_K_A(vector<double> K_A, vector<double> alpha, vector<double> A, vector<double> K_B,
                       vector<double> D_A, vector<double> D_B, vector<double> r, double dr, double M) {
    vector<double> result(alpha.size());
    vector<double> drDalpha(alpha.size());
    vector<double> Dalpha_DB_deriv_sum(alpha.size());
    vector<double> drDB(alpha.size());
    vector<double> D_alpha = spatial_derivative(dr, vec_log(alpha), true, true);
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
                       vector<double> D_A, vector<double> D_B, vector<double> r, double dr, double M) {
    vector<double> result(alpha.size());
    vector<double> dr_ln_psi_vec = dr_ln_psi(r, M);
    vector<double> ddr_ln_psi_vec = ddr_ln_psi(r, M);
    vector<double> D_B_deriv = spatial_derivative(dr, D_B, true, true);
    vector<double> psi_vec = psi(r, M);
    vector<double> D_alpha = spatial_derivative(dr, vec_log(alpha), true, true);


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


vector<double> ddt_alpha(vector<double> alpha, vector<double> K_A, vector<double> K_B) {
    vector<double> result(alpha.size());

    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -2*alpha[i]*(K_A[i]+2*K_B[i]);
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


vector<double> calc_app_horizon(vector<double> A, vector<double> B, vector<double> K_B, vector<double> r, double dr, double M){
    vector<double> result(r.size());
    vector<double> psi_vec = psi(r, M);
    vector<double> B_non_tilde(r.size());

    for (int i = 0; i < r.size(); ++i) {
        B_non_tilde[i] = B[i]*pow(psi_vec[i], 4);
    }

    vector<double> B_deriv = spatial_derivative(dr, B_non_tilde, false, false);

    for (int i = 0; i < r.size(); ++i) {
        result[i] = 1/(sqrt(A[i]*pow(psi_vec[i], 4)))*(2/r[i]+B_deriv[i]/(B[i]*pow(psi_vec[i], 4)))-2*K_B[i];
    }

    return result;
}


double linearInterpolation(const vector<double>& x, const vector<double>& y, double x_0) {
    // Check if arrays have the same size
    if (x.size() != y.size()) {
        cerr << "Error: x and y arrays must have the same size." << endl;
        return 0.0; // Return a default value or handle the error as needed
    }

    // Find the indices of the two points bounding x_0
    size_t i = 0;
    while (i < x.size() && x[i] < x_0) {
        i++;
    }

    // Perform linear interpolation
    if (i == 0) {
        cerr << "Error: x_0 is smaller than the smallest x value." << endl;
        return 0.0; // Return a default value or handle the error as needed
    } else if (i == x.size()) {
        cerr << "Error: x_0 is larger than the largest x value." << endl;
        return 0.0; // Return a default value or handle the error as needed
    } else {
        double x1 = x[i - 1];
        double x2 = x[i];
        double y1 = y[i - 1];
        double y2 = y[i];
        double slope = (y2 - y1) / (x2 - x1);
        double interpolatedY = y1 + slope * (x_0 - x1);
        return interpolatedY;
    }
}


double calc_surf_area(double r_h, vector<double> B, vector<double> r, double M){
    vector<double> psi_vec = psi(r, M);
    vector<double> sh_vec(B.size());

    for (int i = 0; i < sh_vec.size(); ++i) {
        sh_vec[i] = B[i] * pow(psi_vec[i], 4);
    }

    double B_at_rh = linearInterpolation(r, sh_vec, r_h);

    return 4*M_PI*B_at_rh*r_h*r_h;
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
    vector<double> w1_D_A = ddt_D_AB(D_A, alpha, K_A, dr);
    vector<double> w1_D_B = ddt_D_AB(D_B, alpha, K_B, dr);
    vector<double> w1_K_A = ddt_K_A(K_A, alpha, A, K_B, D_A, D_B, r, dr, M);
    vector<double> w1_K_B = ddt_K_B(K_B, alpha, A, B, K_A, D_A, D_B, r, dr, M);
    vector<double> w1_alpha = ddt_alpha(alpha, K_A, K_B);

    vector<double> w2_A = ddt_AB(addv(A, MVS(w1_A, 0.5 * dt)), addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(K_A, MVS(w1_K_A, 0.5 * dt)), r, M);
    vector<double> w2_B = ddt_AB(addv(B, MVS(w1_A, 0.5 * dt)), addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(K_B, MVS(w1_K_B, 0.5 * dt)), r, M);
    vector<double> w2_D_A = ddt_D_AB(addv(D_A, MVS(w1_D_A, 0.5 * dt)), addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(K_A, MVS(w1_K_A, 0.5 * dt)), dr);
    vector<double> w2_D_B = ddt_D_AB(addv(D_B, MVS(w1_D_B, 0.5 * dt)), addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(K_B, MVS(w1_K_B, 0.5 * dt)), dr);
    vector<double> w2_K_A = ddt_K_A(addv(K_A, MVS(w1_K_A, 0.5 * dt)), addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(A, MVS(w1_A, 0.5 * dt)), addv(K_B, MVS(w1_K_B, 0.5 * dt)), addv(D_A, MVS(w1_D_A, 0.5 * dt)), addv(D_B, MVS(w1_D_B, 0.5 * dt)), r, dr, M);
    vector<double> w2_K_B = ddt_K_B(addv(K_B, MVS(w1_K_B, 0.5 * dt)), addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(A, MVS(w1_A, 0.5 * dt)), addv(B, MVS(w1_B, 0.5 * dt)), addv(K_A, MVS(w1_K_A, 0.5 * dt)), addv(D_A, MVS(w1_D_A, 0.5 * dt)),  addv(D_B, MVS(w1_D_B, 0.5 * dt)), r, dr, M);
    vector<double> w2_alpha = ddt_alpha(addv(alpha, MVS(w1_alpha, 0.5 * dt)), addv(K_A, MVS(w1_K_A, 0.5 * dt)), addv(K_B, MVS(w1_K_B, 0.5 * dt)));

    vector<double> w3_A = ddt_AB(addv(A, MVS(w2_A, 0.5 * dt)), addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(K_A, MVS(w2_K_A, 0.5 * dt)), r, M);
    vector<double> w3_B = ddt_AB(addv(B, MVS(w2_A, 0.5 * dt)), addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(K_B, MVS(w2_K_B, 0.5 * dt)), r, M);
    vector<double> w3_D_A = ddt_D_AB(addv(D_A, MVS(w2_D_A, 0.5 * dt)), addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(K_A, MVS(w2_K_A, 0.5 * dt)), dr);
    vector<double> w3_D_B = ddt_D_AB(addv(D_B, MVS(w2_D_B, 0.5 * dt)), addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(K_B, MVS(w2_K_B, 0.5 * dt)), dr);
    vector<double> w3_K_A = ddt_K_A(addv(K_A, MVS(w2_K_A, 0.5 * dt)), addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(A, MVS(w2_A, 0.5 * dt)), addv(K_B, MVS(w2_K_B, 0.5 * dt)), addv(D_A, MVS(w2_D_A, 0.5 * dt)), addv(D_B, MVS(w2_D_B, 0.5 * dt)), r, dr, M);
    vector<double> w3_K_B = ddt_K_B(addv(K_B, MVS(w2_K_B, 0.5 * dt)), addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(A, MVS(w2_A, 0.5 * dt)), addv(B, MVS(w2_B, 0.5 * dt)), addv(K_A, MVS(w2_K_A, 0.5 * dt)), addv(D_A, MVS(w2_D_A, 0.5 * dt)),  addv(D_B, MVS(w2_D_B, 0.5 * dt)), r, dr, M);
    vector<double> w3_alpha = ddt_alpha(addv(alpha, MVS(w2_alpha, 0.5 * dt)), addv(K_A, MVS(w2_K_A, 0.5 * dt)), addv(K_B, MVS(w2_K_B, 0.5 * dt)));

    vector<double> w4_A = ddt_AB(addv(A, MVS(w3_A, dt)), addv(alpha, MVS(w3_alpha, dt)), addv(K_A, MVS(w3_K_A, dt)), r, M);
    vector<double> w4_B = ddt_AB(addv(B, MVS(w3_A, dt)), addv(alpha, MVS(w3_alpha, dt)), addv(K_B, MVS(w3_K_B, dt)), r, M);
    vector<double> w4_D_A = ddt_D_AB(addv(D_A, MVS(w3_D_A, dt)), addv(alpha, MVS(w3_alpha, dt)), addv(K_A, MVS(w3_K_A, dt)), dr);
    vector<double> w4_D_B = ddt_D_AB(addv(D_B, MVS(w3_D_B, dt)), addv(alpha, MVS(w3_alpha, dt)), addv(K_B, MVS(w3_K_B, dt)), dr);
    vector<double> w4_K_A = ddt_K_A(addv(K_A, MVS(w3_K_A, dt)), addv(alpha, MVS(w3_alpha, dt)), addv(A, MVS(w3_A, dt)), addv(K_B, MVS(w3_K_B, dt)), addv(D_A, MVS(w3_D_A, dt)), addv(D_B, MVS(w3_D_B, dt)), r, dr, M);
    vector<double> w4_K_B = ddt_K_B(addv(K_B, MVS(w3_K_B, dt)), addv(alpha, MVS(w3_alpha, dt)), addv(A, MVS(w3_A, dt)), addv(B, MVS(w3_B, dt)), addv(K_A, MVS(w3_K_A, dt)), addv(D_A, MVS(w3_D_A, dt)),  addv(D_B, MVS(w3_D_B, dt)), r, dr, M);
    vector<double> w4_alpha = ddt_alpha(addv(alpha, MVS(w3_alpha, dt)), addv(K_A, MVS(w3_K_A, dt)), addv(K_B, MVS(w3_K_B,  dt)));

    vector<vector<double>> ndata(11, vector<double>(r.size()));

    ndata[0] = addv(A, MVS(addv(addv(addv(w1_A, MVS(w2_A, 2)), MVS(w3_A, 2)), w4_A), dt / 6));
    ndata[1] = addv(B, MVS(addv(addv(addv(w1_B, MVS(w2_B, 2)), MVS(w3_B, 2)), w4_B), dt / 6));
    ndata[2] = addv(D_A, MVS(addv(addv(addv(w1_D_A, MVS(w2_D_A, 2)), MVS(w3_D_A, 2)), w4_D_A), dt / 6));
    ndata[3] = addv(D_B, MVS(addv(addv(addv(w1_D_B, MVS(w2_D_B, 2)), MVS(w3_D_B, 2)), w4_D_B), dt / 6));
    ndata[4] = addv(K_A, MVS(addv(addv(addv(w1_K_A, MVS(w2_K_A, 2)), MVS(w3_K_A, 2)), w4_K_A), dt / 6));
    ndata[5] = addv(K_B, MVS(addv(addv(addv(w1_K_B, MVS(w2_K_B, 2)), MVS(w3_K_B, 2)), w4_K_B), dt / 6));
    ndata[6] = addv(alpha, MVS(addv(addv(addv(w1_alpha, MVS(w2_alpha, 2)), MVS(w3_alpha, 2)), w4_alpha), dt / 6));;
    ndata[7] = spatial_derivative(dr, vec_log(ndata[6]), true, false);
    ndata[8] = r;
    ndata[9] = calc_app_horizon(ndata[0], ndata[1], ndata[5], r, dr, M);

    return ndata;
}


void savetofile(string fname, vector<double> timestepline){
    ofstream of;

    // opening file using ofstream
    of.open(fname, ios::app);
    if (!of){
        ofstream of(fname);
        for (int i = 0; i < timestepline.size()-1; ++i) {
            if (i % 5 == 0) {
                of << fixed << setprecision(8) << timestepline[i] << ", ";
            }
        }
        of << timestepline[timestepline.size()-1] << "\n";
        of.close();
    }
    else {
        for (int i = 0; i < timestepline.size()-1; ++i) {
            if (i % 5 == 0) {
                of << fixed << setprecision(8) << timestepline[i] << ", ";
            }
        }
        of << timestepline[timestepline.size()-1] << "\n";
        of.close();
    }
}



double findZeroPoint(const vector<double>& xData, const vector<double>& yData) {
    if (xData.empty() || yData.empty() || xData.size() != yData.size()) {
        cerr << "Invalid input data!" << endl;
        return 0.0;
    }

    size_t dataSize = xData.size();

    // Find the first data point where x > 0.10 and y changes sign
    size_t firstIndex = 0;
    size_t secondIndex = 0;
    for (size_t i = 1; i < dataSize; ++i) {
        if (xData[i] > 0.10) {
            if ((yData[i - 1] > 0 && yData[i] < 0) || (yData[i - 1] < 0 && yData[i] > 0)) {
                firstIndex = i - 1;
                secondIndex = i;
                break;
            }
        }
    }

    if (secondIndex == 0) {
        cerr << "No intersection point found where x > 0.10!" << endl;
        return 0.0;
    }

    // Perform linear interpolation
    double xFirst = xData[firstIndex];
    double xSecond = xData[secondIndex];
    double yFirst = yData[firstIndex];
    double ySecond = yData[secondIndex];

    double slope = (ySecond - yFirst) / (xSecond - xFirst);
    double xIntersection = xFirst - (yFirst / slope);

    return xIntersection;
}


void timeevolution(const vector<double>& r, double M, double dt, double N_T, double dr){

    vector<vector<double>> tardata(11, vector<double>(r.size()));

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
    tardata[9] = calc_app_horizon(A, B, K_B, r, dr, M);
    savetofile("A.csv", tardata[0]);
    savetofile("B.csv", tardata[1]);
    savetofile("D_A.csv", tardata[2]);
    savetofile("D_B.csv", tardata[3]);
    savetofile("K_A.csv", tardata[4]);
    savetofile("K_B.csv", tardata[5]);
    savetofile("alpha.csv", tardata[6]);
    savetofile("D_alpha.csv", tardata[7]);
    savetofile("r.csv", tardata[8]);
    savetofile("r_BH_apparent.csv", tardata[9]);

    vector<double> r_bh_scalar;
    vector<double> S_bh_scalar;


    for (int i = 0; i < N_T; ++i) {
        tardata = RK4(tardata, dt, M, dr);
        r_bh_scalar.push_back(findZeroPoint(r, tardata[9]));
        S_bh_scalar.push_back(calc_surf_area(r_bh_scalar[r_bh_scalar.size()-1], tardata[1], tardata[8], M));

        if (i % 25 == 0) {
            savetofile("A.csv", tardata[0]);
            savetofile("B.csv", tardata[1]);
            savetofile("D_A.csv", tardata[2]);
            savetofile("D_B.csv", tardata[3]);
            savetofile("K_A.csv", tardata[4]);
            savetofile("K_B.csv", tardata[5]);
            savetofile("alpha.csv", tardata[6]);
            savetofile("D_alpha.csv", tardata[7]);
            savetofile("r.csv", tardata[8]);
            savetofile("r_BH_apparent.csv", tardata[9]);
        }
    }

    savetofile("r_BH_scalar.csv", r_bh_scalar);
    savetofile("S_BH_apparent.csv", S_bh_scalar);
}

using namespace chrono;

int main() {
    auto start = high_resolution_clock::now();
    filesystem::path cwd = filesystem::current_path();

    vector<string> filenames = {
            "A.csv",
            "B.csv",
            "D_A.csv",
            "D_B.csv",
            "K_A.csv",
            "K_B.csv",
            "alpha.csv",
            "D_alpha.csv",
            "r.csv",
            "r_BH_apparent.csv",
            "S_BH_apparent.csv",
            "r_BH_scalar.csv",
    };

    for (const string& filename : filenames) {
        filesystem::path filepath = cwd / filename;
        if (filesystem::exists(filepath)) {
            if (filesystem::remove(filepath)) {
                cout << "File deleted successfully: " << filepath << endl;
            } else {
                cerr << "Error deleting the file: " << filepath << endl;
                return 1;
            }
        } else {
            cout << "File does not exist: " << filepath << endl;
        }
    }

    double dr = 0.01;
    vector<double> r = generate_linspace(10000, dr);

    timeevolution(r, 1., 0.005, 2000, dr);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    auto elapsed = static_cast<float>(duration.count());
    cout << "Code took " << elapsed/1000000 << " seconds to execute" << endl;

    return 1;
}