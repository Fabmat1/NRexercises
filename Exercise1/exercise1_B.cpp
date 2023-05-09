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


vector<double> euler_timestepper (int n_timesteps, int gridsize){

    vector<double> x(gridsize);
    x = generate_linspace(gridsize);
    double gridfraction = 1./gridsize;

    vector<double> y(gridsize);
    y = u(x);

    vector<double> y_prime(gridsize);
    y_prime= central_diff(y, gridfraction);

    vector<double> u_step(x.size());
    vector<double> RMSE(n_timesteps);

    vector<double> timestepped_x(x.size());
    timestepped_x = x;
    vector<double> analytic_sol(x.size());
    double this_step;
    double timefraction = 1./n_timesteps;

    for (int i = 0; i < n_timesteps; ++i) {
        for (int j = 0; j < x.size(); ++j) {
            timestepped_x[j] += timefraction;
        }
        analytic_sol = u(timestepped_x);

        for (int k = 0; k < x.size(); ++k) {
            this_step = y[k] + timefraction*y_prime[k];
            u_step[k] = this_step;
            RMSE[i] += pow(this_step-analytic_sol[k],2);
        }
        RMSE[i] = sqrt(gridfraction*RMSE[i]);
        y_prime = central_diff(u_step, gridfraction);
        y = u_step;
    }

    return RMSE;
}


void example_RMSE(int n_timesteps, int n_gridpoints){
    vector<double> RMSE = euler_timestepper(n_timesteps, n_gridpoints);

    cout << "Writing output..." << endl;
    vector<double> timestep_vector(RMSE.size());
    # pragma omp parallel for
    for (int i = 0; i < RMSE.size(); ++i) {
        timestep_vector[i] = (1./n_timesteps)*i;
    }
    save_file(timestep_vector, RMSE, "RMSE.txt");
}


void save_matrix(vector<vector<double>> matrix,  const std::string &filepath){
    ofstream outputFile(filepath);

    for (auto & i : matrix) {
        for (int j = 0; j < i.size(); ++j) {
            outputFile << i[j] << "\t";
            if(j == i.size()-1){
                outputFile << "\n";
            }
        }
    }

    outputFile.close();
    std::cout << "Saved output." << std::endl;
}


int main() {
    cout << "Generating data for B.2" << endl;
    int gridsize = 1000;
    double gridfraction = 1./gridsize;
    vector<double> x = generate_linspace(gridsize);
    vector<double> y = u(x);
    vector<double> y_prime = central_diff(y, gridfraction);

    save_file(x, y_prime, "b2_data.txt");

    cout << "Done with B.2" << endl;

    example_RMSE(50000, 1000);

    cout << "Probing Euler method timestepper for 10 to 100 gridpoints and 10 to 100 timesteps" << endl;


    vector<vector<double>> RMSE_matrix(90, vector<double>(90));
    for (int i = 0; i < 90; ++i) {
        for (int j = 0; j < 90; ++j) {
            cout << i << j;
            vector<double> RMSE(i+10);
            RMSE = euler_timestepper(i+10, j+10);
            double matrix_entry = RMSE.back();
            RMSE_matrix[i][j] = matrix_entry;
        }
    }

//    save_matrix(RMSE_matrix, "RMSEspace.txt");

    cout << "All done!" << endl;
    return 1;
}