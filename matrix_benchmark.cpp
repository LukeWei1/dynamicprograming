// Benchmark. Build with `make bench`.
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "matrix.h"

template <class F> double medianMs(F f, int reps) {
    std::vector<double> t;
    for (int i = 0; i < reps; ++i) {
        auto s = std::chrono::steady_clock::now();
        auto r = f();
        auto e = std::chrono::steady_clock::now();
        if (r->size() == 0) std::cerr << "";  // keep result alive
        t.push_back(std::chrono::duration<double, std::milli>(e - s).count());
    }
    std::sort(t.begin(), t.end());
    return t[t.size() / 2];
}

int main() {
    std::mt19937_64 rng(42);
    const std::vector<std::size_t> cutoffs = {16, 32, 64, 128, 256};
    std::cout << "Median runtime in ms (3 runs). Strassen columns show the cutoff size.\n";
    std::cout << std::setw(6) << "n" << std::setw(12) << "conv";
    for (auto c : cutoffs) std::cout << std::setw(10) << ("S/" + std::to_string(c));
    std::cout << std::setw(10) << "pureS" << "\n";
    for (std::size_t n : {64, 128, 256, 512, 1024, 2048}) {
        Matrix a = Matrix::random(n, rng), b = Matrix::random(n, rng);
        std::cout << std::setw(6) << n << std::fixed << std::setprecision(1)
                  << std::setw(12) << medianMs([&] { return a.multiply(b); }, 3);
        for (auto c : cutoffs)
            std::cout << std::setw(10) << medianMs([&] { return a.strassenMultiply(b, c); }, 3);
        if (n <= 512) std::cout << std::setw(10) << medianMs([&] { return a.strassenMultiply(b, 1); }, 3);
        std::cout << std::endl;
    }
}
