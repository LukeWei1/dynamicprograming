// matrix.h -- square matrix with conventional and Strassen multiplication.
#pragma once
#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <vector>

/// A square n x n matrix of doubles stored in row-major order.
class Matrix {
public:
    /// Size at or below which Strassen falls back to conventional multiplication.
    static constexpr std::size_t DEFAULT_CUTOFF = 64;

    /// Creates an n x n zero matrix (n must be > 0).
    explicit Matrix(std::size_t n);

    /// Builds a matrix from rows; the number of rows must equal each row's length.
    static Matrix fromRows(const std::vector<std::vector<double>>& rows);
    /// The n x n identity matrix.
    static Matrix identity(std::size_t n);
    /// An n x n matrix with entries uniformly random in [-1, 1).
    static Matrix random(std::size_t n, std::mt19937_64& rng);
    /// Reassembles four equal-sized quadrants (inverse of split()).
    static Matrix join(const Matrix& c11, const Matrix& c12, const Matrix& c21, const Matrix& c22);

    std::size_t size() const { return n_; }

    /// Bounds-checked element access (throws std::out_of_range).
    double get(std::size_t r, std::size_t c) const;
    /// Bounds-checked element assignment (throws std::out_of_range).
    void set(std::size_t r, std::size_t c, double v);
    double& operator()(std::size_t r, std::size_t c);
    double operator()(std::size_t r, std::size_t c) const;

    /// Splits into quadrants {A11, A12, A21, A22}. Size must be even.
    std::vector<Matrix> split() const;
    /// Copy padded with zeros to m x m (m >= n).
    Matrix padded(std::size_t m) const;
    /// Top-left k x k block (1 <= k <= n).
    Matrix cropped(std::size_t k) const;

    Matrix operator+(const Matrix& o) const;
    Matrix operator-(const Matrix& o) const;
    /// Operator form of multiply(): `auto c = a * b;` (empty if sizes differ).
    std::optional<Matrix> operator*(const Matrix& o) const;

    /// Conventional O(n^3) multiplication. Empty optional if dimensions are incompatible.
    std::optional<Matrix> multiply(const Matrix& o) const;

    /// Strassen multiplication. Non-power-of-two sizes are zero-padded and the result cropped.
    /// Sub-problems of size <= cutoff use conventional multiplication (hybrid); cutoff = 1 is pure Strassen.
    /// Empty optional if dimensions are incompatible.
    std::optional<Matrix> strassenMultiply(const Matrix& o, std::size_t cutoff = DEFAULT_CUTOFF) const;

    bool operator==(const Matrix& o) const;
    /// True if all entries agree within a relative tolerance eps.
    bool approxEquals(const Matrix& o, double eps = 1e-9) const;
    std::string toString() const;

private:
    std::size_t n_;
    std::vector<double> data_;

    void check(std::size_t r, std::size_t c) const;
    void requireSame(const Matrix& o) const;
    Matrix multiplyConventional(const Matrix& o) const;
    Matrix strassen(const Matrix& o, std::size_t cutoff) const;  // n must be a power of two
};
