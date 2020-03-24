#include <iostream>

#include "matrix.h"

template<unsigned N>
Matrix<N, 1>
solve(const Matrix<N, N> &lhs, const Matrix<N, 1> &rhs)
{
    return lhs.inverse() * rhs;
}

template<unsigned N>
Matrix<N, N>
identity()
{
    return Matrix<N, N>::identity();
}

template<unsigned N>
bool
is_invertible(Matrix <N, N> m)
{
    double det = m.determinant();
    return -0.1 > det || det > 0.1;
}

int
main(void) {
    Matrix<2U, 2U> matrix_a {
        2., 1.,
        4., 2. };

    auto id = identity<2>();
    bool has_inverse = is_invertible(id);
    std::cout << "identity is inversible: " << (has_inverse ? "yes\n" : "no\n");

    Matrix<2U, 2U> matrix_b = matrix_a + (5 ^ id);
    Matrix<2, 1> res = solve(matrix_b,  { 6., 10. });
    std::cout << "RESULT IS:\n" << res;

    return 0;
    (void) has_inverse;
}