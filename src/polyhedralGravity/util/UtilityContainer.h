#pragma once

#include <array>
#include <set>
#include <numeric>
#include <utility>
#include <algorithm>
#include <functional>
#include <cmath>
#include <string>
#include <iostream>

namespace polyhedralGravity::util {

    /**
     * Alias for a two-dimensional array with size M and N. M is the major size!
     */
    template<typename T, size_t M, size_t N>
    using Matrix = std::array<std::array<T, N>, M>;

    /**
     * Applies a binary function to elements of two containers piece by piece. The objects must
     * be iterable and should have the same size!
     * @tparam Container an iterable object like an array or vector
     * @tparam BinOp a binary function to apply
     * @param lhs the first container
     * @param rhs the second container
     * @param binOp a binary function like +, -, *, /
     * @return a container containing the result
     */
    template<typename Container, typename BinOp>
    Container applyBinaryFunction(const Container &lhs, const Container &rhs, BinOp binOp) {
        Container ret = lhs;
        std::transform(std::begin(lhs), std::end(lhs), std::begin(rhs), std::begin(ret), binOp);
        return ret;
    }

    /**
     * Applies a binary function to elements of one container piece by piece. The objects must
     * be iterable. The resulting container consists of the containers' object after the application
     * of the binary function with the scalar as a parameter.
     * @tparam Container an iterable object like an array or vector
     * @tparam Scalar a scalar to use on each element
     * @tparam BinOp a binary function to apply
     * @param lhs the first container
     * @param scalar a scalar to use on each element
     * @param binOp a binary function like +, -, *, /
     * @return a container containing the result
     */
    template<typename Container, typename Scalar, typename BinOp>
    Container applyBinaryFunction(const Container &lhs, const Scalar &scalar, BinOp binOp) {
        Container ret = lhs;
        std::transform(std::begin(lhs), std::end(lhs), std::begin(ret), [&binOp, &scalar](const Scalar &element) {
            return binOp(element, scalar);
        });
        return ret;
    }

    /**
     * Applies the Operation Minus to two Containers piece by piece.
     * @code {1, 2, 3} - {1, 1, 1} = {0, 1, 2} @endcode
    * @tparam Container an iterable object like an array or vector
     * @param lhs minuend
     * @param rhs subtrahend
     * @return the difference
     */
    template<typename Container>
    Container operator-(const Container &lhs, const Container &rhs) {
        return applyBinaryFunction(lhs, rhs, std::minus<>());
    }

    /**
    * Applies the Operation Plus to two Containers piece by piece.
    * @code {1, 2, 3} + {1, 1, 1} = {2, 3, 4} @endcode
    * @tparam Container an iterable object like an array or vector
    * @param lhs addend
    * @param rhs addend
    * @return the sum
    */
    template<typename Container>
    Container operator+(const Container &lhs, const Container &rhs) {
        return applyBinaryFunction(lhs, rhs, std::plus<>());
    }

    /**
    * Applies the Operation * to two Containers piece by piece.
    * @code {1, 2, 3} * {2, 2, 2} = {2, 4, 6} @endcode
    * @tparam Container an iterable object like an array or vector
    * @param lhs multiplicand
    * @param rhs multiplicand
    * @return the product
    */
    template<typename Container>
    Container operator*(const Container &lhs, const Container &rhs) {
        return applyBinaryFunction(lhs, rhs, std::multiplies<>());
    }

    /**
    * Applies the Operation / to two Containers piece by piece.
    * @code {1, 2, 3} / {1, 2, 3} = {1, 1, 1} @endcode
    * @tparam Container an iterable object like an array or vector
    * @param lhs multiplicand
    * @param rhs multiplicand
    * @return the product
    */
    template<typename Container>
    Container operator/(const Container &lhs, const Container &rhs) {
        return applyBinaryFunction(lhs, rhs, std::divides<>());
    }

    /**
    * Applies the Operation + to a Container and a Scalar.
    * @code {1, 2, 3} + 2 = {3, 4, 5} @endcode
    * @tparam Container an iterable object like an array or vector
    * @tparam Scalar a scalar to use on each element
    * @param lhs addend
    * @param scalar addend
    * @return a Container
    */
    template<typename Container, typename Scalar>
    Container operator+(const Container &lhs, const Scalar &scalar) {
        return applyBinaryFunction(lhs, scalar, std::plus<>());
    }

    /**
    * Applies the Operation to a Container and a Scalar.
    * @code {1, 2, 3} - 2 = {-1, 0, 1} @endcode
    * @tparam Container an iterable object like an array or vector
    * @tparam Scalar a scalar to use on each element
    * @param lhs minuend
    * @param scalar subtrahend
    * @return a Container
    */
    template<typename Container, typename Scalar>
    Container operator-(const Container &lhs, const Scalar &scalar) {
        return applyBinaryFunction(lhs, scalar, std::minus<>());
    }

    /**
    * Applies the Operation to a Container and a Scalar.
    * @code {1, 2, 3} * 2 = {2, 4, 6} @endcode
    * @tparam Container an iterable object like an array or vector
    * @tparam Scalar a scalar to use on each element
    * @param lhs multiplicand
    * @param scalar multiplicand
    * @return a Container
    */
    template<typename Container, typename Scalar>
    Container operator*(const Container &lhs, const Scalar &scalar) {
        return applyBinaryFunction(lhs, scalar, std::multiplies<>());
    }

    /**
     * Applies the Operation / to a Container and a Scalar.
     * @code {2, 4, 6} / 2 = {1, 2, 3} @endcode
     * @tparam Container an iterable object like an array or vector
     * @tparam Scalar a scalar to use on each element
     * @param lhs the dividend
     * @param scalar the divisor
     * @return a Container
     */
    template<typename Container, typename Scalar>
    Container operator/(const Container &lhs, const Scalar &scalar) {
        return applyBinaryFunction(lhs, scalar, std::divides<>());
    }

    /**
     * Applies the Euclidean norm/ L2-norm to a Container (e.g., a vector)
     * @tparam Container must be iterable
     * @param container e.g. a vector
     * @return a double containing the L2 norm
     */
    template<typename Container>
    double euclideanNorm(const Container &container) {
        return std::sqrt(std::inner_product(std::begin(container), std::end(container), std::begin(container), 0.0));
    }

    /**
     * Computes the absolute value for each value in the given container
     * @tparam Container an iterable container, containing numerical values
     * @param container the container
     * @return a container with the modified values
     */
    template<typename Container>
    Container abs(const Container &container) {
        Container ret = container;
        std::transform(std::begin(container), std::end(container), std::begin(ret),
                       [](const auto &element) { return std::abs(element); });
        return ret;
    }

    /**
     * Computes the determinant with the Sarrus rule for a 3x3 matrix.
     * Notice that for square matrices @f$det(A) = det(A^T)@f$.
     * @tparam T a numerical value
     * @param matrix the 3x3 matrix
     * @return the determinant
     */
    template<typename T>
    T det(const Matrix<T, 3, 3> &matrix) {
        return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] +
               matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0] -
               matrix[0][0] * matrix[1][2] * matrix[2][1] - matrix[0][1] * matrix[1][0] * matrix[2][2];
    }

    /**
     * Computes the transposed of a mxn matrix.
     * @tparam T the type of the matrix elements
     * @tparam M the row number
     * @tparam N the column number
     * @param matrix the matrix to transpose
     * @return the transposed
     */
    template<typename T, size_t M, size_t N>
    Matrix<T, M, N> transpose(const Matrix<T, M, N> &matrix) {
        Matrix<T, N, M> transposed;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                transposed[i][j] = matrix[j][i];
            }
        }
        return transposed;
    }

    /**
    * Returns the cross-product of two cartesian vectors.
    * @tparam T a number
    * @param lhs left vector
    * @param rhs right vector
    * @return cross product
    */
    template<typename T>
    std::array<T, 3> cross(const std::array<T, 3> &lhs, const std::array<T, 3> &rhs) {
        std::array<T, 3> result{};
        result[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
        result[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
        result[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
        return result;
    }

    /**
     * Calculates the normal N as (first * second) / |first * second| with * being the cross product and first, second
     * as cartesian vectors.
     * @tparam T numerical type
     * @param first first cartesian vector
     * @param second second cartesian vector
     * @return the normal (normed)
     */
    template<typename T>
    std::array<T, 3> normal(const std::array<T, 3> &first, const std::array<T, 3> &second) {
        const std::array<T, 3> crossProduct = cross(first, second);
        const double norm = euclideanNorm(crossProduct);
        return crossProduct / norm;
    }

    /**
    * Returns the dot product of two cartesian vectors.
    * @tparam T a number
    * @param lhs left vector
    * @param rhs right vector
    * @return dot product
    */
    template<typename T>
    T dot(const std::array<T, 3> &lhs, const std::array<T, 3> &rhs) {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
    }

    /**
     * Implements the signum function with a certain EPSILON to absorb rounding errors.
     * @tparam T a numerical (floating point) value
     * @param val the value itself
     * @param cutoffEpsilon the cut-off radius around zero to return 0
     * @return -1, 0, 1 depending on the sign and the given EPSILON
     */
    template<typename T>
    int sgn(T val, double cutoffEpsilon) {
        return val < -cutoffEpsilon ? -1 : val > cutoffEpsilon ? 1 : 0;
    }

    /**
     * Concatenates two std::array of different sizes to one array.
     * @tparam T the shared type of the arrays
     * @tparam M the size of the first container
     * @tparam N  the size of the second container
     * @param first the first array
     * @param second the second array
     * @return a new array of size M+N with type T
     */
    template<typename T, size_t M, size_t N>
    std::array<T, M + N> concat(const std::array<T, M> &first, const std::array<T, N> &second) {
        std::array<T, M + N> result{};
        size_t index = 0;
        for (const auto &el: first) {
            result[index++] = el;
        }
        for (const auto &el: second) {
            result[index++] = el;
        }
        return result;
    }

    /**
     * Calculates the surface area of a triangle consisting of three cartesian vertices.
     * @tparam T numerical type
     * @return surface area
     */
    template<typename T>
    T surfaceArea(const Matrix<T, 3, 3> &triangle) {
        return 0.5 * euclideanNorm(cross(triangle[1] - triangle[0], triangle[2] - triangle[0]));
    }

    /**
     * Calculates the magnitude between two values and return true if the magnitude
     * between the exponents in greater than 17.
     * @tparam T numerical type
     * @param first first number
     * @param second second number
     * @return true if the difference is too huge, so that floating point absorption will happen
     */
    template<typename T>
    bool isCriticalDifference(const T &first, const T &second) {
        // 50 is the (log2) exponent of the floating point (1 / 1e-15)
        constexpr int maxExponentDifference = 50;
        int x, y;
        std::frexp(first, &x);
        std::frexp(second, &y);
        return std::abs(x - y) > maxExponentDifference;
    }

    /**
    * A utility struct that creates an overload set out of all the function objects it inherits from.
    * It allows a lambda function to be used with std::visit in a variant.
    * The lambda function needs to be able to handle all types contained in the variant,
    * and this struct allows it to do so by inheriting the function-call operator from each type.
    * @tparam Ts A template parameter pack representing all types for which the struct should be able to handle.
    *
    * @ref https://en.cppreference.com/w/cpp/utility/variant/visit
    */
    template<class... Ts>
    struct overloaded : Ts... {
        using Ts::operator()...;
    };

    /**
    * This function template provides deduction guides for the overloaded struct.
    * It deduces the template arguments for overloaded based on its constructor arguments
    * thus saving you from explicitly mentioning them while instantiation.
    * @tparam Ts A template parameter pack representing all types that will be passed to the overloaded struct.
    * @param Ts Variable numbers of parameters to pass to the overloaded struct.
    * @return This doesn't return a value, but rather assists in the creation of an overloaded object.
    *          The type of the object will be overloaded<Ts...>.
    *
    * @ref https://en.cppreference.com/w/cpp/utility/variant/visit
    */
    template<class... Ts>
    overloaded(Ts...) -> overloaded<Ts...>;

    namespace detail {

        /**
         * Helper method for the tuple operator+, which expands the tuple into a parameter pack.
         * @tparam Ts types of the tuple
         * @tparam Is indices of the tuple
         * @param t1 first tuple
         * @param t2 second tuple
         * @return a new tuple with the added values
         */
        template<typename... Ts, size_t... Is>
        auto tuple_add(const std::tuple<Ts...> &t1, const std::tuple<Ts...> &t2, std::index_sequence<Is...>) {
            return std::make_tuple(std::get<Is>(t1) + std::get<Is>(t2)...);
        }

    }

    /**
     * Adds the contents of two tuples of the same size and types with the operator +.
     * @tparam Ts types of the tuples
     * @param t1 first tuple
     * @param t2 second tuple
     * @return a new tuple with the added values
     */
    template<typename... Ts>
    auto operator+(const std::tuple<Ts...> &t1, const std::tuple<Ts...> &t2) {
        return detail::tuple_add(t1, t2, std::index_sequence_for<Ts...>{});
    }

    /**
     * Operator << for an array of any size.
     * @tparam T type of the array, must have an << operator overload
     * @tparam N size of the array
     * @param os the ostream
     * @param array the array itself
     * @return ostream
     */
    template<typename T, size_t N>
    std::ostream &operator<<(std::ostream &os, const std::array<T, N> &array) {
        os.operator<<('[');
        os.operator<<(' ');
        std::for_each(array.cbegin(), array.cend(), [&os](const auto& arg) {
            os << arg << ' ';
        });
        os.operator<<(']');
        return os;
    }

    /**
     * Operator << for a set.
     * @tparam T type of the set, must have an << operator overload
     * @param os the ostream
     * @param set the set itself
     * @return ostream
     */
    template<typename T>
    std::ostream &operator<<(std::ostream &os, const std::set<T> &set) {
        os.operator<<('[');
        os.operator<<(' ');
        std::for_each(set.cbegin(), set.cend(), [&os](const auto& arg) {
            os << arg << ' ';
        });
        os.operator<<(']');
        return os;
    }

    template<typename T>
    struct is_stdarray : std::false_type {
    };

    template<typename T, std::size_t N>
    struct is_stdarray<std::array<T, N>> : std::true_type {
    };

}