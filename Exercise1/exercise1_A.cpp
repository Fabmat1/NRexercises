#include "exercise1_A.h"
#include "vector"
#include "omp.h"
#include <immintrin.h>
#include "chrono"
#include "iostream"
#include <fstream>

using namespace std::chrono;
using namespace std;

// Function that generates a "logarithmic" space such that x_i = x_{i-1}+h and x_{i+1} = x_i + rh for all x
vector<double> generate_logspace(int length, double h, double r, double start = 0.0) {
    vector<double> result(length);
    result[0] = start;
    result[1] = start + h;

    for (int i = 2; i < length; i++) {
        result[i] = result[i - 1] + (result[i - 1] - result[i - 2]) * r;
    }

    return result;
}


vector<double> cubic_function(vector<double> x, double a = 1., double b = 0., double c = 0., double d = 0.) {
    const int n = x.size();
    vector<double> y(n);

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] = x[i] * x[i] * x[i] * a + x[i] * x[i] * b + c * x[i] + d; // ax^3+bx^2+cx+d
    }
    return y;
}


pair<vector<double>, vector<double>> three_point_stencil(const vector<double> &x_arr, const vector<double> &y_arr) {
    const int n = y_arr.size() - 2;
    vector<double> dydx(n);
    std::vector<double> truncated_x_arr(n);

    #pragma omp parallel for
    for (int i = 1; i < n; i++) {
        dydx[i] = (y_arr[i + 1] - y_arr[i - 1]) /
                  (x_arr[i + 1] - x_arr[i - 1])-((x_arr[i+1]-x_arr[i])*(x_arr[i+1]-x_arr[i])-(x_arr[i]-x_arr[i-1])*(x_arr[i]-x_arr[i-1]))/((x_arr[i + 1] - x_arr[i - 1])); // (f(x_(i+1))-f_(x_(i-1)))/(x_(i+1)-x_(i-1))
        truncated_x_arr[i] = x_arr[i];
    }
    return {truncated_x_arr, dydx};
}


void save_file(vector<double> vector1, vector<double> vector2, const string &filepath = "../output.txt") {
    ofstream outputFile(filepath);

    for (int i = 0; i < vector1.size(); i++) {
        outputFile << vector1[i] << "\t" << vector2[i] << std::endl;
    }

    outputFile.close();
    std::cout << "Saved output." << std::endl;
}


int main() {
    auto start = high_resolution_clock::now();

    vector<double> logspace = generate_logspace(1000, 1e-10, 1.1, 0);
    vector<double> y = cubic_function(logspace, 0, 1, 0, 0);
    save_file(logspace, y, "../pure_function.txt");

    auto derivative = three_point_stencil(logspace, y);
    save_file(derivative.first, derivative.second);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    auto elapsed = static_cast<float>(duration.count());
    std::cout << "Code took " << elapsed / 1000 << "ms to execute." << std::endl;
    std::cout << "All done! " << std::endl;
    return 1;
}