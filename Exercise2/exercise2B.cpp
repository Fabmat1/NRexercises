#include "vector"
#include "omp.h"
#include <immintrin.h>
#include "chrono"
#include "iostream"
#include <fstream>
#include "cmath"
#include <functional>
#include <utility>

using namespace std::chrono;
using namespace std;


vector<double> spatial_derivative(double dx, const vector<double> &y, bool neg_bound) {
    int s = y.size();
    int factor;
    if (neg_bound) {
        factor = -1;
    } else {
        factor = 1;
    }

    vector<double> result(s);
    for (int i = 2; i < y.size() - 2; ++i) {
        result[i] = (1. / 12 * y[i - 2] - 2./3 * y[i - 1] + 2./3 * y[i + 1] - 1. / 12 * y[i + 2]) / (dx);
    }
    result[0] = (1. / 12 * y[2] * factor - 2./3 * y[1] * factor + 2./3 * y[1] - 1. / 12 * y[2]) / (dx);
    result[1] = (1. / 12 * y[1] * factor - 2./3 * y[0] + 2./3 * y[2] - 1. / 12 * y[3]) / (dx);
    result[s - 1] = (1. / 12 * y[s - 3] - 2./3 * y[s - 2] + 2./3 * y[s - 2] - 1. / 12 * y[s - 3]) / (dx);
    result[s - 2] = (1. / 12 * y[s - 4] - 2./3 * y[s - 3] + 2./3 * y[s - 1] - 1. / 12 * y[s - 1]) / (dx);
    return result;
}


vector<double> ddt_AB(vector<double> AB, vector<double> alpha, vector<double> K_AB) {
    vector<double> result(alpha.size());
    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -2 * alpha[i] * AB[i] * K_AB[i];
    }
    return result;
}

function<vector<double>(vector<double>)> RK_friendly_dt_AB(const function<vector<double>(vector<double>, vector<double>,
                        vector<double>)>& A, const vector<double>& alpha, const vector<double>& K_AB) {
    return [=](vector<double> arg) {
        return A(arg, alpha, K_AB);
    };
}


vector<double>
ddt_D_AB(const vector<double>& D_AB, vector<double> alpha, vector<double> K_AB, vector<double> D_alpha, double dr) {
    vector<double> result(alpha.size());
    vector<double> K_AB_deriv = spatial_derivative(dr, K_AB, false);
    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -2 * alpha[i] * (K_AB[i] * D_alpha[i] + K_AB_deriv[i]);
    }
    return result;
}

function<vector<double>(vector<double>)> RK_friendly_dt_D_AB(
        const function<vector<double>(vector<double>, vector<double>,vector<double>,vector<double>,double)>& A,
        const vector<double>& alpha, const vector<double>& K_AB, vector<double> D_alpha, double dr) {
    return [=](vector<double> arg) {
        return A(arg, alpha, K_AB, D_alpha, dr);
    };
}


vector<double> ddt_K_A(vector<double> K_A, vector<double> alpha, vector<double> A, vector<double> K_B,
                       vector<double> D_A, vector<double> D_B, vector<double> D_alpha, vector<double> r, double dr) {
    vector<double> result(alpha.size());
    vector<double> Dalpha_DB_sum(D_alpha.size());
    for (int i = 0; i < D_alpha.size(); ++i) {
        Dalpha_DB_sum[i] = D_alpha[i] + D_B[i];
    }
    vector<double> Dalpha_DB_sum_deriv(alpha.size());
    Dalpha_DB_sum_deriv = spatial_derivative(dr, Dalpha_DB_sum, true);

    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -alpha[i] / A[i] * (Dalpha_DB_sum_deriv[i] + D_alpha[i] * D_alpha[i] - (D_alpha[i] * D_B[i]) / 2
                                        + (D_B[i] * D_B[i]) / 2 - (D_A[i] * D_B[i]) / 2 -
                                        A[i] * K_A[i] * (K_A[i] + 2 * K_B[i]) - 1 / r[i] * (D_A[i] - 2 * D_B[i]));
    }
    return result;
}

function<vector<double>(vector<double>)> RK_friendly_dt_K_A(
        const function<vector<double>(vector<double>, vector<double>, vector<double>, vector<double>,
                                      vector<double>, vector<double>, vector<double>, vector<double>, double)>& F,
        vector<double> alpha, vector<double> A, vector<double> K_B,
        vector<double> D_A, vector<double> D_B, vector<double> D_alpha, vector<double> r, double dr) {
    return [=](vector<double> arg) {
        return F(arg, alpha, A, K_B, D_A, D_B, D_alpha, r, dr);
    };
}


vector<double> ddt_K_B(vector<double> K_B,vector<double> alpha, vector<double> A, vector<double> B, vector<double> K_A,
                       vector<double> D_A, vector<double> D_B, vector<double> D_alpha, vector<double> r, double dr) {
    vector<double> result(alpha.size());
    vector<double> D_B_deriv = spatial_derivative(dr, D_B, true);

    for (int i = 0; i < alpha.size(); ++i) {
        result[i] = -alpha[i] / (2 * A[i]) * (D_B_deriv[i] + D_alpha[i] * D_B[i] + D_B[i] * D_B[i] - (D_A[i] * D_B[i]) / 2
            - 1 / r[i] * (D_A[i] - 2 * D_alpha[i] - 4 * D_B[i]) -
            (2 * (A[i] - B[i])) / (r[i] * r[i] * B[i]))
            + alpha[i] * K_B[i] * (K_A[i] + 2 * K_B[i]);
    }
    return result;
}


function<vector<double>(vector<double>)> RK_friendly_dt_K_B(
        const function<vector<double>(vector<double>,vector<double>, vector<double>, vector<double>, vector<double>,
                                      vector<double>, vector<double>, vector<double>, vector<double>, double)>& F,
        vector<double> alpha, vector<double> A, vector<double> B, vector<double> K_A,
        vector<double> D_A, vector<double> D_B, vector<double> D_alpha, vector<double> r, double dr) {
    return [=](vector<double> arg) {
        return F(arg, alpha, A, B, K_A, D_A, D_B, D_alpha, r, dr);
    };
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


vector<double> RK4(vector<double> data, const function<vector<double>(vector<double>)> derivative, double dt) {
    vector<double> w1 = derivative(data);
    vector<double> w2 = derivative(addv(data, MVS(w1, 0.5 * dt)));
    vector<double> w3 = derivative(addv(data, MVS(w2, 0.5 * dt)));
    vector<double> w4 = derivative(addv(data, MVS(w2, dt)));
    vector<double> ndata =
            addv(data, MVS(addv(addv(addv(w1, MVS(w2, 2)), MVS(w3, 2)), w4), dt / 6));
    return ndata;
}


void timeevolution(vector<double> r, double M_BH, double dr, int N_timestep, double dt, vector<double> initial_alpha,
                   vector<double> initial_K_A, vector<double> initial_K_B) {
    vector<double> phi(r.size());
    vector<double> initial_D_A(r.size());
    vector<double> initial_D_B(r.size());
    vector<double> initial_A(r.size());
    vector<double> initial_B(r.size());
    vector<double> initial_D_alpha(r.size());
    for (int i = 0; i < phi.size(); ++i) {
        phi[i] = 1 + M_BH / (2 * r[i]);
    }
    initial_A = phi;
    initial_B = phi;

    vector<vector<double>> initial_D_arrays;
    for (const vector<double> &Arr: {initial_A, initial_B, initial_alpha}) {
        vector<double> temp(r.size());
        for (int i = 0; i < r.size(); ++i) {
            temp[i] = log(Arr[i]);
        }
        temp = spatial_derivative(dr, temp, false);

        initial_D_arrays.push_back(temp);
    }
    initial_D_A = initial_D_arrays[0];
    initial_D_B = initial_D_arrays[1];
    initial_D_alpha = initial_D_arrays[2];

    vector<double> step_D_A(r.size());
    vector<double> step_D_B(r.size());
    vector<double> step_A(r.size());
    vector<double> step_B(r.size());
    vector<double> step_K_A(r.size());
    vector<double> step_K_B(r.size());

    step_A = RK4(initial_A, RK_friendly_dt_AB(ddt_AB, initial_alpha, initial_K_A), dt);
    step_B = RK4(initial_B, RK_friendly_dt_AB(ddt_AB, initial_alpha, initial_K_A), dt);
    step_D_A = RK4(initial_D_A, RK_friendly_dt_D_AB(ddt_D_AB, initial_alpha, initial_K_A, initial_D_alpha, dr), dt);
    step_D_B = RK4(initial_D_B, RK_friendly_dt_D_AB(ddt_D_AB, initial_alpha, initial_K_B, initial_D_alpha, dr), dt);
    step_K_A = RK4(initial_K_A, RK_friendly_dt_K_A(ddt_K_A, initial_alpha, initial_A, initial_K_B, initial_D_A, initial_D_B, initial_D_alpha, r, dr), dt);
    step_K_B = RK4(initial_K_B, RK_friendly_dt_K_B(ddt_K_B, initial_alpha, initial_A, initial_B, initial_K_A, initial_D_A, initial_D_B, initial_D_alpha, r, dr), dt);

    for (int i = 1; i < N_timestep; ++i) {
        step_A = RK4(step_A, RK_friendly_dt_AB(ddt_AB, initial_alpha, step_K_A), dt);
        step_B = RK4(step_B, RK_friendly_dt_AB(ddt_AB, initial_alpha, step_K_A), dt);
        step_D_A = RK4(step_D_A, RK_friendly_dt_D_AB(ddt_D_AB, initial_alpha, step_K_A, initial_D_alpha, dr), dt);
        step_D_B = RK4(step_D_B, RK_friendly_dt_D_AB(ddt_D_AB, initial_alpha, step_K_B, initial_D_alpha, dr), dt);
        step_K_A = RK4(step_K_A, RK_friendly_dt_K_A(ddt_K_A, initial_alpha, step_A, step_K_B, step_D_A, step_D_B, initial_D_alpha, r, dr), dt);
        step_K_B = RK4(step_K_B, RK_friendly_dt_K_B(ddt_K_B, initial_alpha, step_A, step_B, step_K_A, step_D_A, step_D_B, initial_D_alpha, r, dr), dt);
    }

    for (double i : step_A) {
        cout << i << endl;
    }

}


int main() {
    vector<double> r = generate_linspace(1000, 0.01);
    vector<double> rprime = spatial_derivative(0.01, r, false);


    vector<double> alpha(r.size(), 1);
    vector<double> K_A(r.size(), 1);
    vector<double> K_B(r.size(), 1);

    timeevolution(r, 100000., 0.01, 2000, 0.005, alpha, K_A, K_B);

    return 1;
}