#include "vector"
#include "omp.h"
#include <immintrin.h>
#include "chrono"
#include "iostream"
#include <fstream>

using namespace std::chrono;
using namespace std;


// Function that generates a "logarithmic" space such that x_i = x_{i-1}+h and x_{i+1} = x_i + rh for all i
vector<double> generate_logspace(int length, double h0, double r, double start = 0.0) {
    vector<double> result(length);
    result[0] = start;
    result[1] = start + h0;

    // This cannot be parallelized. Sad.
    for (int i = 2; i < length; i++) {
        result[i] = result[i - 1] + (result[i - 1] - result[i - 2]) * r;
    }

    return result;
}


vector<double> polynomial(vector<double> x, double a = 1., double b = 0., double c = 0., double d = 0.) {
    const int n = x.size();
    vector<double> y(n);

    // Use multiprocessing for fast loop execution
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] = x[i] * x[i] * x[i] * a + x[i] * x[i] * b + c * x[i] + d; // ax^3+bx^2+cx+d
    }
    return y;
}


// Implementation of the non-uniform three-point stencil formula
pair<vector<double>, vector<double>> non_uniform_three_point_stencil(const vector<double> &x_arr, const vector<double> &y_arr, double r) {
    const int n = y_arr.size() - 2;
    vector<double> dydx(n);
    std::vector<double> truncated_x_arr(n);

    // Use multiprocessing for fast loop execution
    #pragma omp parallel for
    for (int i = 1; i < n; i++) {
        // Refer to PDF to compare formulas
        dydx[i] = ((y_arr[i]- y_arr[i-1])*r*r + y_arr[i + 1]-y_arr[i]) /((x_arr[i] - x_arr[i - 1])*r*(1+r));
        truncated_x_arr[i] = x_arr[i];
    }
    return {truncated_x_arr, dydx};
}


// Implementation of the uniform second derivative three-point stencil formula
pair<vector<double>, vector<double>> sd_three_point_stencil(const vector<double> &x_arr, const vector<double> &y_arr, double r) {
    const int n = y_arr.size() - 2;
    vector<double> ddydxx(n);
    std::vector<double> truncated_x_arr(n);

    // Use multiprocessing for fast loop execution
    #pragma omp parallel for
    for (int i = 1; i < n+1; i++) {
        // Refer to PDF to compare formulas
        ddydxx[i-1] = 2*((y_arr[i-1]+ (1/r)*y_arr[i+1])-(1+1/r)*y_arr[i]) /((1+r)*(x_arr[i] - x_arr[i - 1])*(x_arr[i] - x_arr[i - 1]));
        truncated_x_arr[i-1] = x_arr[i];
    }

    return {truncated_x_arr, ddydxx};
}


void save_file(vector<double> vector1, vector<double> vector2, const std::string &filepath) {
    ofstream outputFile(filepath);

    for (int i = 0; i < vector1.size(); i++) {
        outputFile << vector1[i] << "\t" << vector2[i] << "\n";
    }

    outputFile.close();
    std::cout << "Saved output." << std::endl;
}


int main() {
    auto start = high_resolution_clock::now();

    std::cout << "Working on task A.1" << std::endl;
    double r = 1.1;

    vector<double> logspace = generate_logspace(1000, 1e-10, r, 0);
    vector<double> y = polynomial(logspace, 0, 1, 0, 0);

    auto derivative = non_uniform_three_point_stencil(logspace, y, r);
    // Save for plotting with matplotlib
    save_file(derivative.first, derivative.second, "first_derivative.txt");

    std::cout << "Finished working on task A.1, starting task A.2" << std::endl;

    vector<double> z = polynomial(logspace, 0, 1, 0, 0);


    auto second_derivative =sd_three_point_stencil(logspace, z, r);
    save_file(second_derivative.first, second_derivative.second, "second_derivative.txt");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    auto elapsed = static_cast<float>(duration.count());
    std::cout << "Code took " << elapsed / 1000 << "ms to execute." << std::endl;
    std::cout << "All done! " << std::endl;
    return 1;
}