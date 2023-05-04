#include "vector"
#include "omp.h"
#include <immintrin.h>
#include "chrono"
#include "iostream"
#include <fstream>
#include "cmath"

using namespace std::chrono;
using namespace std;

// Function that generates a linear space in the interval (0, 1) such that x_i = i/N
vector<double> generate_linspace(int length) {
    double fl_length = static_cast< double >(length);
    vector<double> result(length);

    // Use AVX2 instructions for efficiency
    __m256d vec_l = _mm256_set_pd(fl_length, fl_length, fl_length, fl_length);

    // Vectorized and parallelized loop
    #pragma omp parallel for
    for (int i = 0; i < length; i += 4) {
        __m256d vec_i = _mm256_set_pd(i+3, i+2, i+1, i);

        __m256d vec_result = _mm256_div_pd(vec_i, vec_l);
        _mm256_storeu_pd(&result[i], vec_result);
    }

    if (length % 4 != 0){
        for (int i = length-4; i < length; ++i) {
            result[i] = i/fl_length;
        }
    }

    return result;
}



// Implementation of the finite central difference approximation
vector<double>  central_diff(const vector<double> &y_arr, double dx) {
    const int n = y_arr.size();
    vector<double> dydx(n);

    // Use multiprocessing for fast loop execution
    #pragma omp parallel for
    for (int i = 1; i < n-1; i++) {
        // Refer to PDF to compare formulas
        dydx[i] = (y_arr[i+1]-y_arr[i-1]) /(2*dx);
    }

    dydx[0] =  (y_arr[1]-y_arr[y_arr.size()-1]) /(2*dx);
    dydx[y_arr.size()-1] =  (y_arr[0]-y_arr[y_arr.size()-2]) /(2*dx);

    return dydx;
}


vector<double> u(const vector<double>& x){
    vector<double> result(x.size());

    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i) {
        result[i] = exp(-2*cos(2*numbers::pi*x[i]));
    }

    return result;
}

void save_file(vector<double> vector1, vector<double> vector2, const std::string &filepath) {
    ofstream outputFile(filepath);
    if (vector1.size() != vector2.size()){
        cout << "Vector size mismatch" << endl;
        return;
    }

    for (int i = 0; i < vector1.size(); i++) {
        outputFile << vector1[i] << "\t" << vector2[i] << "\n";
    }

    outputFile.close();
    std::cout << "Saved output." << std::endl;
}

int main() {
    cout << "Generating data for B.2" << endl;
    int gridsize = 1000;
    auto gridfraction = static_cast<double>(gridsize);
    gridfraction = 1/gridfraction;
    vector<double> x = generate_linspace(gridsize);
    vector<double> y = u(x);
    vector<double> y_prime = central_diff(y, gridfraction);

    save_file(x, y_prime, "b2_data.txt");

    cout << "Done with B.2" << endl;

    vector<double> u_step(x.size());
    int n_timesteps = 100000;
    auto timefraction = static_cast<double>(n_timesteps);
    timefraction = 1/timefraction;

    vector<double> RMSE(n_timesteps);
    vector<double> timestepped_x(x.size());
    vector<double> analytic_sol(x.size());
    double this_step;

    cout << "Using Euler method timestepper with " << n_timesteps << " steps" << endl;
    timestepped_x = x;
    for (int i = 0; i < n_timesteps; ++i) {
        # pragma omp parallel for
        for (int j = 0; j < x.size(); ++j) {
            timestepped_x[j] += timefraction;
        }
        analytic_sol = u(timestepped_x);

        # pragma omp parallel for
        for (int k = 0; k < x.size(); ++k) {
            this_step = y[k] + timefraction*y_prime[k];
            u_step[k] = this_step;
            RMSE[i] += pow(this_step-analytic_sol[k],2);
        }
        RMSE[i] = sqrt(gridfraction*RMSE[i]);
        y_prime = central_diff(u_step, gridfraction);
        y = u_step;
    }

    cout << "Writing output..." << endl;
    vector<double> timestep_vector(RMSE.size());
    # pragma omp parallel for
    for (int i = 0; i < RMSE.size(); ++i) {
        timestep_vector[i] = timefraction*i;
    }
    save_file(timestep_vector, RMSE, "RSME.txt");

    cout << "All done!" << endl;
    return 1;
}