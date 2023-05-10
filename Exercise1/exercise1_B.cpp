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
        __m256d vec_i = _mm256_set_pd(i + 3, i + 2, i + 1, i);

        __m256d vec_result = _mm256_div_pd(vec_i, vec_l);
        _mm256_storeu_pd(&result[i], vec_result);
    }

    if (length % 4 != 0) {
        for (int i = length - 4; i < length; ++i) {
            result[i] = i / fl_length;
        }
    }

    return result;
}


// Implementation of the finite central difference approximation
vector<double> central_diff(const vector<double> &y_arr, double dx) {
    const int n = y_arr.size();
    vector<double> dydx(n);

    // Use multiprocessing for fast loop execution
#pragma omp parallel for
    for (int i = 1; i < n - 1; i++) {
        // Refer to PDF to compare formulas
        dydx[i] = (y_arr[i + 1] - y_arr[i - 1]) / (2 * dx);
    }

    dydx[0] = (y_arr[1] - y_arr[y_arr.size() - 1]) / (2 * dx);
    dydx[y_arr.size() - 1] = (y_arr[0] - y_arr[y_arr.size() - 2]) / (2 * dx);

    return dydx;
}


vector<double> u(const vector<double> &x) {
    vector<double> result(x.size());
#pragma omp parallel for
    for (int i = 0; i < x.size(); ++i) {
        result[i] = exp(-2 * cos(2 * M_PI * x[i]));
    }

    return result;
}

double interpolate_y_prime(const vector<double>& x, const vector<double>& y, double p) {
    while (p > 1){
        p -= 1;
    }

    while (p < 0){
        p += 1;
    }

    int n = x.size();

    size_t i = 0;
    while (i < n && x[i] < p) {
        i++;
    }

    // Perform linear interpolation
    double x0 = x[i - 1];
    double x1 = x[i];
    double y0 = y[i - 1];
    double y1 = y[i];

    double slope = (y1 - y0) / (x1 - x0);
    double interpolated_y = y0 + slope * (p - x0);

    return interpolated_y;
}


void save_file(vector<double> vector1, vector<double> vector2, const std::string &filepath) {
    ofstream outputFile(filepath);
    if (vector1.size() != vector2.size()) {
        cout << "Vector size mismatch" << endl;
        return;
    }

    for (int i = 0; i < vector1.size(); i++) {
        outputFile << vector1[i] << "\t" << vector2[i] << "\n";
    }

    outputFile.close();
    std::cout << "Saved output to " << filepath << std::endl;
}


vector<double> euler_timestepper(int n_timesteps, int gridsize) {
    vector<double> RMSE(n_timesteps);
    double RMSE_proxy = 0;
    vector<double> x(gridsize);
    vector<double> y_step(gridsize);
    vector<double> x_aly(gridsize);
    vector<double> y_aly(gridsize);
    vector<double> y = u(x);
    vector<double> y_prime = central_diff(y, 1. / gridsize);

    for (int i = 0; i < x.size(); ++i) {
        x[i] = i * 1. / gridsize;
    }
    y = u(x);

    for (int i = 0; i < n_timesteps; ++i) {
#pragma omp parallel for
        for (int j = 0; j < gridsize; ++j) {
            x_aly[j] = x[j] + i * 1. / n_timesteps;
            y_step[j] = y[j] + 1. / n_timesteps * y_prime[j];
        }

        y_aly = u(x_aly);

        for (int j = 0; j < gridsize; ++j) {
            RMSE_proxy += (y[j] - y_aly[j]) * (y[j] - y_aly[j]);
        }

        y = y_step;
        y_prime = central_diff(y, 1. / gridsize);

        RMSE[i] = sqrt((1. / gridsize) * RMSE_proxy);
        RMSE_proxy = 0;

    }

    return RMSE;
}


void example_RMSE(int n_timesteps, int n_gridpoints) {
    vector<double> RMSE = euler_timestepper(n_timesteps, n_gridpoints);

    cout << "Writing output..." << endl;
    vector<double> timestep_vector(RMSE.size());
# pragma omp parallel for
    for (int i = 0; i < RMSE.size(); ++i) {
        timestep_vector[i] = (1. / n_timesteps) * i;
    }
    save_file(timestep_vector, RMSE, "RMSE.txt");
}


void save_matrix(vector<vector<double>> matrix, const std::string &filepath) {
    ofstream outputFile(filepath);

    for (auto &i: matrix) {
        for (int j = 0; j < i.size(); ++j) {
            if (j == i.size() - 1) {
                outputFile << i[j] << "\n";
            } else {
                outputFile << i[j] << "\t";
            }
        }
    }

    outputFile.close();
    std::cout << "Saved Table to " << filepath << std::endl;
}


double avg(const vector<double> &arr) {
    double mean;

    for (double i: arr) {
        mean += i;
    }
    mean /= arr.size();

    return mean;
}


vector<double> RK4_timestepper(int n_timesteps, int gridsize) {
    vector<double> RMSE(n_timesteps);
    double RMSE_proxy = 0;
    vector<double> x(gridsize);
    vector<double> y_step(gridsize);
    vector<double> x_aly(gridsize);
    vector<double> y_aly(gridsize);
    double w1;
    double w2;
    double w3;
    double w4;

    for (int i = 0; i < x.size(); ++i) {
        x[i] = i * 1. / gridsize;
    }

    vector<double> y = u(x);
    vector<double> y_prime = central_diff(y, 1. / gridsize);

    for (int i = 0; i < n_timesteps; ++i) {
        #pragma omp parallel for
        for (int j = 0; j < gridsize; ++j) {
            x_aly[j] = x[j] + i * 1. / n_timesteps;
            w1 = y_prime[j];
            w2 = interpolate_y_prime(x, y_prime, x[j] + w1/(2*n_timesteps));
            w3 = interpolate_y_prime(x, y_prime, x[j] + w2/(2*n_timesteps));
            w4 = interpolate_y_prime(x, y_prime, x[j] + w3/n_timesteps);
            y_step[j] = y[j] + (w1 + 2 * w2 + 2 * w3 + w4) / (6 * n_timesteps) ;
//            y_step[j] = y[j] + 1. / n_timesteps * w1;
        }

        y_aly = u(x_aly);

        for (int j = 0; j < gridsize; ++j) {
            RMSE_proxy += (y[j] - y_aly[j]) * (y[j] - y_aly[j]);
        }

        y = y_step;
        y_prime = central_diff(y, 1. / gridsize);

        RMSE[i] = sqrt((1. / gridsize) * RMSE_proxy);
        RMSE_proxy = 0;

        if (i == 15000){
            save_file(x, y, "test_RK4.txt");
        }
    }

    return RMSE;
}


void example_RMSE_RK4(int n_timesteps, int n_gridpoints) {
    vector<double> RMSE = RK4_timestepper(n_timesteps, n_gridpoints);

    cout << "Writing output..." << endl;
    vector<double> timestep_vector(RMSE.size());
# pragma omp parallel for
    for (int i = 0; i < RMSE.size(); ++i) {
        timestep_vector[i] = (1. / n_timesteps) * i;
    }
    save_file(timestep_vector, RMSE, "RMSE_RK4.txt");
}


void probe_timestepper(vector<double> stepper (int, int), int N_t, int N_x, const string& saveloc, int stepsize_t=10, int stepsize_x=10){

    cout << "Probing timestepper for 10 to "<< N_x*stepsize_x <<" gridpoints and 10 to "<< N_t*stepsize_t <<" timesteps" << endl;

    vector<vector<double>> RMSE_matrix(N_t, vector<double>(N_x));
    for (int i = 0; i < N_t; i++) {
        for (int j = 0; j < N_x; j++) {
            vector<double> RMSE(i + 10);
            RMSE = stepper(stepsize_t*i + 10, stepsize_x*j + 10);
            RMSE_matrix[i][j] = avg(RMSE);
        }
    }

    save_matrix(RMSE_matrix, saveloc);

    ofstream outputFile(saveloc+"_params.txt");
    outputFile << N_t*stepsize_t << "\t";
    outputFile << N_x*stepsize_x << "\t";

    outputFile.close();
}


int main() {
    cout << "Generating data for B.2" << endl;
    int gridsize = 1000;
    double gridfraction = 1. / gridsize;
    vector<double> x = generate_linspace(gridsize);
    vector<double> y = u(x);
    vector<double> y_prime = central_diff(y, gridfraction);

    save_file(x, y_prime, "b2_data.txt");

    cout << "Done with B.2, starting B.3" << endl;

    example_RMSE(100000, 1000);
    probe_timestepper(euler_timestepper, 20, 20, "RMSEspace.txt", 25, 5);

    cout << "Starting B.4 by redoing tasks from B.3 using RK4" << endl;

    example_RMSE_RK4(100000, 1000);
    probe_timestepper(RK4_timestepper, 20, 20, "RMSEspace_RK4.txt", 25, 5);

    cout << "All done!" << endl;
    return 1;
}