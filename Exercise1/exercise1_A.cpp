#include "exercise1_A.h"
#include "vector"
#include "omp.h"
#include <immintrin.h>
#include "chrono"
#include "iostream"

using namespace std::chrono;
using namespace std;

// Function that generates a "logarithmic" space such that x_i = x_{i-1}+h and x_{i+1} = x_i + rh for all x
vector<double> generate_logspace(int length, double h, double r, double start=0.0){
    vector<double> result(length);
    result[0] = start;
    result[1] = start+h;

    for (int i = 2; i < length; i++){
        result[i] = result[i-1]+(result[i-1]-result[i-2])*r;
    }

    return result;
}


int main(){
    auto start = high_resolution_clock::now();

    vector<double> logspace = generate_logspace(1000, 1e-10, 1.1, 0);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    auto elapsed = static_cast<float>(duration.count());
    std::cout << elapsed/1000 << std::endl;
    return 1;
}