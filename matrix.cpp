// matrix.cpp -- implementation of Matrix (see matrix.h for documentation).
#include "matrix.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

Matrix::Matrix(std::size_t n) : n_(n), data_(n * n, 0.0) {}

//takes nested vectors. Creates Matrix m where number of rows is the size
//if row length is different from count throws error
//inner loop writes data into each index and returns m
Matrix Matrix::fromRows(const std::vector<std::vector<double>>& rows) {
    Matrix m(rows.size());
    for (std::size_t r = 0; r < rows.size(); ++r) {
        if (rows[r].size() != rows.size()) throw std::invalid_argument("Matrix must be square");
        for (std::size_t c = 0; c < rows.size(); ++c) m.data_[r * m.n_ + c] = rows[r][c];
    }
    return m;
}

//creates an identity matrix of size n
Matrix Matrix::identity(std::size_t n) {
    Matrix m(n);
    for (std::size_t i = 0; i < n; ++i) m.data_[i * n + i] = 1.0;
    return m;
}

//dist maps random number generator output between [-1, 1]
//loops by references sets each element to dist
Matrix Matrix::random(std::size_t n, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Matrix m(n);
    for (double& x : m.data_) x = dist(rng);
    return m;
}
//throws error if row or column is too big.
void Matrix::check(std::size_t r, std::size_t c) const {
    if (r >= n_ || c >= n_) throw std::out_of_range("Matrix index out of range");
}
//checks if matricies are the same size
void Matrix::requireSame(const Matrix& o) const {
    if (n_ != o.n_) throw std::invalid_argument("Matrix size mismatch");
}
//get and set validate index then read/write
double Matrix::get(std::size_t r, std::size_t c) const { check(r, c); return data_[r * n_ + c]; }
void Matrix::set(std::size_t r, std::size_t c, double v) { check(r, c); data_[r * n_ + c] = v; }
double& Matrix::operator()(std::size_t r, std::size_t c) { check(r, c); return data_[r * n_ + c]; }
double Matrix::operator()(std::size_t r, std::size_t c) const { return get(r, c); }

//splits matrix up
//odd size can't be cut in half so throws error
std::vector<Matrix> Matrix::split() const {
    if (n_ % 2 != 0) throw std::invalid_argument("Can only split even-sized matrices");
    //h is the quadrant size
    const std::size_t h = n_ / 2;
    //q holds four h×h zero matrices, in the order top-left, top-right, bottom-left, bottom-right.
    std::vector<Matrix> q(4, Matrix(h));
    for (std::size_t r = 0; r < h; ++r) {
        std::copy_n(&data_[r * n_], h, &q[0].data_[r * h]); //A11
        std::copy_n(&data_[r * n_ + h], h, &q[1].data_[r * h]); //A12
        std::copy_n(&data_[(r + h) * n_], h, &q[2].data_[r * h]);//A13
        std::copy_n(&data_[(r + h) * n_ + h], h, &q[3].data_[r * h]);//A14
    }
    return q;
}
//joins matricies

Matrix Matrix::join(const Matrix& c11, const Matrix& c12, const Matrix& c21, const Matrix& c22) {
    const std::size_t h = c11.n_, n = 2 * h;
    if (c12.n_ != h || c21.n_ != h || c22.n_ != h) throw std::invalid_argument("Quadrant size mismatch");
    Matrix res(n);
    for (std::size_t r = 0; r < h; ++r) {
        std::copy_n(&c11.data_[r * h], h, &res.data_[r * n]);
        std::copy_n(&c12.data_[r * h], h, &res.data_[r * n + h]);
        std::copy_n(&c21.data_[r * h], h, &res.data_[(r + h) * n]);
        std::copy_n(&c22.data_[r * h], h, &res.data_[(r + h) * n + h]);
    }
    return res;
}
//Makes a bigger zero matrix and copies each of this matrix’s rows into its top-left corner.
Matrix Matrix::padded(std::size_t m) const {
    if (m < n_) throw std::invalid_argument("Cannot pad to a smaller size");
    Matrix res(m);
    for (std::size_t r = 0; r < n_; ++r) std::copy_n(&data_[r * n_], n_, &res.data_[r * m]);
    return res;
}
//oppisite of padding
Matrix Matrix::cropped(std::size_t k) const {
    if (k == 0 || k > n_) throw std::invalid_argument("Invalid crop size");
    Matrix res(k);
    for (std::size_t r = 0; r < k; ++r) std::copy_n(&data_[r * n_], k, &res.data_[r * k]);
    return res;
}

Matrix Matrix::operator+(const Matrix& o) const {
    requireSame(o);
    Matrix res(n_);
    for (std::size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] + o.data_[i];
    return res;
}

Matrix Matrix::operator-(const Matrix& o) const {
    requireSame(o);
    Matrix res(n_);
    for (std::size_t i = 0; i < data_.size(); ++i) res.data_[i] = data_[i] - o.data_[i];
    return res;
}

std::optional<Matrix> Matrix::operator*(const Matrix& o) const { return multiply(o); }

std::optional<Matrix> Matrix::multiply(const Matrix& o) const {
    if (o.n_ != n_) return std::nullopt;
    return multiplyConventional(o);
}

// conventional matrix multiplication
Matrix Matrix::multiplyConventional(const Matrix& o) const {
    Matrix res(n_);
    for (std::size_t i = 0; i < n_; ++i)
        for (std::size_t k = 0; k < n_; ++k) {
            const double a = data_[i * n_ + k];
            const double* brow = &o.data_[k * n_];
            double* crow = &res.data_[i * n_];
            for (std::size_t j = 0; j < n_; ++j) crow[j] += a * brow[j];
        }
    return res;
}
//multiplication using strassen method
std::optional<Matrix> Matrix::strassenMultiply(const Matrix& o, std::size_t cutoff) const {
    if (o.n_ != n_) return std::nullopt;
    if (cutoff < 1) throw std::invalid_argument("cutoff must be >= 1");
    std::size_t m = 1;
    //Finds the smallest power of two that is as greater than or equal to size
    while (m < n_) m *= 2;
    if (m == n_) return strassen(o, cutoff);
    //pads matricies if needed
    return padded(m).strassen(o.padded(m), cutoff).cropped(n_);
}

// Seven products follow the Wikipedia formulation of Strassen's algorithm.
Matrix Matrix::strassen(const Matrix& o, std::size_t cutoff) const {
    if (n_ <= cutoff || n_ == 1) return multiplyConventional(o);
    auto a = split();
    auto b = o.split();
    Matrix m1 = (a[0] + a[3]).strassen(b[0] + b[3], cutoff);
    Matrix m2 = (a[2] + a[3]).strassen(b[0], cutoff);
    Matrix m3 = a[0].strassen(b[1] - b[3], cutoff);
    Matrix m4 = a[3].strassen(b[2] - b[0], cutoff);
    Matrix m5 = (a[0] + a[1]).strassen(b[3], cutoff);
    Matrix m6 = (a[2] - a[0]).strassen(b[0] + b[1], cutoff);
    Matrix m7 = (a[1] - a[3]).strassen(b[2] + b[3], cutoff);
    return join(m1 + m4 - m5 + m7, m3 + m5, m2 + m4, m1 - m2 + m3 + m6);
}

bool Matrix::operator==(const Matrix& o) const { return n_ == o.n_ && data_ == o.data_; }

bool Matrix::approxEquals(const Matrix& o, double eps) const {
    if (n_ != o.n_) return false;
    for (std::size_t i = 0; i < data_.size(); ++i) {
        double scale = std::max({1.0, std::fabs(data_[i]), std::fabs(o.data_[i])});
        if (std::fabs(data_[i] - o.data_[i]) > eps * scale) return false;
    }
    return true;
}

std::string Matrix::toString() const {
    std::string s;
    for (std::size_t r = 0; r < n_; ++r)
        for (std::size_t c = 0; c < n_; ++c) s += std::to_string(data_[r * n_ + c]) + (c + 1 < n_ ? " " : "\n");
    return s;
}
