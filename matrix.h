#include <algorithm>
#include <iostream>
#include <limits>

#define UNUSED(x) ((void)x)

template <bool b>
struct Do_not_compile_when;

template <>
struct Do_not_compile_when<false>
{
    Do_not_compile_when() {}
};

template <>
struct Do_not_compile_when<true>
{
    Do_not_compile_when() = delete;
};

#define FAIL_COMPILE_IF(COND) do {(void)Do_not_compile_when<(COND)>();} while(0)
#define FAIL_COMPILE FAIL_COMPILE_IF(false)

namespace matrix_tools {
// Helper for the matrix constructor
    template<unsigned N, unsigned M, unsigned CurrOffset, typename ...Args>
    struct Helper {
        static void helper(double (&array)[N][M], Args ...args) {
            UNUSED(array);
            UNUSED(sizeof...(args));
            FAIL_COMPILE;
        }
    };

    template<unsigned N, unsigned M, unsigned CurrOffset, typename ...Args>
    struct Helper<N, M, CurrOffset, double, Args...> {
        static void helper(double (&array)[N][M], double x, Args ...args) {
            FAIL_COMPILE_IF(CurrOffset >= N * M);
            array[CurrOffset / M][CurrOffset % M] = x;
            Helper<N, M, CurrOffset + 1U, Args...>::helper(array, args...);
        }
    };

    template<unsigned N, unsigned M>
    struct Helper<N, M, N * M> {
        static void helper(double (&array)[N][M]) {
            UNUSED(array);
            return;
        }
    };

    template<unsigned N, unsigned M, typename...Args>
    void fill_array(double (&array)[N][M], double first, Args ...args)
    {
        Helper<N, M, 0, double, Args...>::helper(array, first, args...);
    }
}

template <unsigned N, unsigned M,
          template <unsigned A, unsigned B> class Parent>
    class Matrix_base;

template <
template <unsigned A, unsigned B> class Parent,
    unsigned I, unsigned J, unsigned K>
    Parent<I, K>
    operator *(
        const Matrix_base<I, J, Parent> &m1,
        const Matrix_base<J, K, Parent> &m2);

template <unsigned N, unsigned M,
          template <unsigned A, unsigned B> class Parent>
class Matrix_base {

public:
    double content[N][M] = {};

public:
    Matrix_base() {}

    template<typename ...Args>
    Matrix_base(double first, Args ...args) {
        matrix_tools::fill_array<N, M>(content, first, args...);
    }

private:

    //@ requires x < N && y < M;
    double get(unsigned x, unsigned y) {
        return content[x][y];
    }

    //@ requires orig < N && dest < M;
    void
    add_line(unsigned orig, unsigned dest, double scalar)
    {
        for (unsigned j = 0; j < M; ++j) {
            content[dest][j] = content[dest][j] + scalar * content[orig][j];
        }
    }

    //@ requires lineno < N;
    void
    scalar_line_multiplication(unsigned lineno, double scalar)
    {
        for (unsigned j = 0; j < M; ++j) {
            content[lineno][j] = scalar * content[lineno][j];
        }
    }

public:

    friend
    std::ostream &
    operator<<(std::ostream& out, const Matrix_base &mat)
    {
        out << "<" << N << ", " << M << ">\n";
        for (unsigned i = 0; i < N; i++) {
            for(unsigned j = 0; j < M; j++) {
                out << mat.content[i][j] << " ";
            }
            out << "\n";
        }
        return out;
    }

    friend
    Parent<N, M>
    operator^(double scalar, const Parent<N, M> &mat)
    {
        Parent<N, M> res;
        for (unsigned i = 0; i < N; i++) {
            for(unsigned j = 0; j < M; j++) {
                res.content[i][j] =
                    scalar * mat.content[i][j];
            }
        }
        return res;
    }

    friend
    Parent<N, M>
    operator +(
        const Matrix_base<N, M, Parent> &m1,
        const Matrix_base<N, M, Parent> &m2)
    {
        Parent<N, M> res;
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                res.Matrix_base<N, M, Parent>::content[i][j] =
                    m1.content[i][j] + m2.content[i][j];
            }
        }
        return res;
    }

    template <
        template <unsigned A, unsigned B> class Res_Parent,
        unsigned I, unsigned J, unsigned K>
    friend
    Res_Parent<I, K>
    operator *(
        const Matrix_base<I, J, Parent> &m1,
        const Matrix_base<J, K, Parent> &m2);

    // Template to SFINAE it when N != M
    template<bool b = true>
    double
    trace() const
    {
        FAIL_COMPILE_IF(N != M);
        double res = 0;
        for(unsigned k = 0; k < N; ++k)
            res+= content[k][k];
        return res;
    }

    // Template to SFINAE it when N != M
    template<bool b = true>
    double
    determinant() const
    {
        FAIL_COMPILE_IF(N != M);
        auto diag = [&](unsigned shift) -> double
        {
            double res = 1;
            for (unsigned i = 0; i < N; ++i) {
                res = res * content[i % N][(i + shift) % N];
            }
            return res;
        };

        double res = 0;
        for (unsigned i = 0; i < N; ++i) {
            if (i%2 == 0)
                res += diag(i);
            else
                res -= diag(i);
        }
        return res;
    }

    // Template to SFINAE it when N != M
    template<bool b = true>
    Parent<N, N>
    inverse() const
    {
        // ultra-naive gaussian elimination.

        FAIL_COMPILE_IF(N != M);
        Parent<N, N> res = Parent<N, N>::identity();
        Parent<N, N> current = *this;

        // Diagonalize current.
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                if (i == j)
                    continue;
                double scalar = -current.get(j, i) / current.get(i, i);
                current.add_line(i, j, scalar);
                res.add_line(i, j, scalar);
            }
        }

        // Normalize diagonal. Do it only by acting on res, since we do not care
        // anymore for res.
        for (unsigned i = 0; i < N; ++i) {
            double scalar = 1 / current.Matrix_base::content[i][i];
            // current.scalar_line_multiplication(i, scalar);
            res.scalar_line_multiplication(i, scalar);
        }

        return res;
    }

    double
    norm_1() const
    {
        double min = std::numeric_limits<double>::lowest();
        for (unsigned i = 0; i < N; ++i) {
            double current = 0;
            for(unsigned j = 0; j < M; ++j) {
                current += content[i][j];
            }
            min = std::min(min, current);
        }
        return min;
    }
};

template <
template <unsigned A, unsigned B> class Parent,
    unsigned I, unsigned J, unsigned K>
    Parent<I, K>
    operator *(
        const Matrix_base<I, J, Parent> &m1,
        const Matrix_base<J, K, Parent> &m2)
{
    auto cross_product =
        [&m1, &m2] (unsigned i, unsigned j) -> double
        {
            double sum = 0;
            for (unsigned k = 0; k < J; ++k) {
                sum += m1.content[i][k] * m2.content[k][j];
            }
            return sum;
        };

    Parent<I, K> res;
    for (unsigned i = 0; i < I; i++) {
        for(unsigned j = 0; j < K; j++) {
            res.content[i][j] = cross_product(i, j);
        }
    }
    return res;
}


template<unsigned N, unsigned M>
class Matrix: public Matrix_base<N, M, Matrix>
{
public:
    template<typename ...Args>
    Matrix(Args ...args): Matrix_base<N, M, Matrix>(args...) {}

    Matrix<N, 1>
    operator [](unsigned n) {
        FAIL_COMPILE_IF(N == 1);
        Matrix <N, 1> res;
        for (unsigned i=0; i < N; ++i) {
            res.content[i][0] = Matrix_base<N, M, Matrix>::content[i][n];
        }
        return res;
    }

    static
    Matrix<N, N>
    identity() {
        Matrix<N, N> res;
        for (unsigned i = 0; i < N; ++i) {
            res.Matrix_base<N, N, Matrix>::content[i][i] = 1;
        }

        return res;
    }
};

template<>
class Matrix<1, 1>: public Matrix_base<1, 1, Matrix>
{
public:
    Matrix(double d): Matrix_base<1, 1, Matrix>(d) {}

    operator double() { return Matrix_base<1, 1, Matrix>::content[0][0]; }
};

template<unsigned N>
class Matrix<N, 1>: public Matrix_base<N, 1, Matrix>
{
public:
    template<typename ...Args>
    Matrix(Args ...args): Matrix_base<N, 1, Matrix>(args...) {}

    Matrix<1, 1>
    operator [](unsigned n) {
        return { Matrix_base<N, 1, Matrix>::content[n][0] };
    }

    operator double() { return Matrix_base<N, 1, Matrix>::content[0][0]; }
};