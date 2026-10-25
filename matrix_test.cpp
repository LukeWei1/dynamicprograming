// Unit tests for Matrix. Build with `make test`.
#include <iostream>
#include "matrix.h"

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; std::cerr << "FAIL line " << __LINE__ << ": " #cond "\n"; } } while (0)

template <class F> bool throws(F f) { try { f(); } catch (...) { return true; } return false; }

int main() {
    // get / set / bounds
    Matrix m(3);
    CHECK(m.get(1, 2) == 0.0);
    m.set(1, 2, 7.5);
    CHECK(m(1, 2) == 7.5);
    CHECK(throws([&] { m.get(3, 0); }));
    

    // split / join
    std::mt19937_64 rng(1);
    Matrix r8 = Matrix::random(8, rng);
    auto q = r8.split();
    CHECK(q[0].size() == 4);
    CHECK(q[0].get(0, 0) == r8.get(0, 0));
    CHECK(q[1].get(0, 0) == r8.get(0, 4));
    CHECK(q[2].get(0, 0) == r8.get(4, 0));
    CHECK(q[3].get(3, 3) == r8.get(7, 7));
    CHECK(Matrix::join(q[0], q[1], q[2], q[3]) == r8);
    CHECK(throws([] { Matrix(3).split(); }));

    // pad / crop
    Matrix r5 = Matrix::random(5, rng);
    CHECK(r5.padded(8).get(7, 7) == 0.0);
    CHECK(r5.padded(8).cropped(5) == r5);

    // + and -
    Matrix a = Matrix::fromRows({{1, 2}, {3, 4}});
    Matrix b = Matrix::fromRows({{5, 6}, {7, 8}});
    CHECK(a + b == Matrix::fromRows({{6, 8}, {10, 12}}));
    CHECK(b - a == Matrix::fromRows({{4, 4}, {4, 4}}));

    // known products
    Matrix expected = Matrix::fromRows({{19, 22}, {43, 50}});
    CHECK(*a.multiply(b) == expected);
    CHECK(*(a * b) == expected);
    CHECK(*a.strassenMultiply(b, 1) == expected);

    // identity
    Matrix r10 = Matrix::random(10, rng);
    CHECK(r10.approxEquals(*r10.multiply(Matrix::identity(10))));
    CHECK(r10.approxEquals(*r10.strassenMultiply(Matrix::identity(10), 1)));

    // incompatible sizes
    CHECK(!Matrix(2).multiply(Matrix(4)).has_value());
    CHECK(!Matrix(2).strassenMultiply(Matrix(4)).has_value());

    // Strassen == conventional on many sizes / cutoffs (incl. non powers of two)
    for (std::size_t n : {1, 2, 3, 4, 5, 7, 8, 16, 17, 31, 33, 64, 100}) {
        Matrix x = Matrix::random(n, rng), y = Matrix::random(n, rng);
        Matrix ref = *x.multiply(y);
        for (std::size_t cutoff : {1, 2, 4, 16}) {
            CHECK(ref.approxEquals(*x.strassenMultiply(y, cutoff), 1e-8));
        }
    }

    if (failures == 0) std::cout << "All Matrix tests passed\n";
    return failures == 0 ? 0 : 1;
}
