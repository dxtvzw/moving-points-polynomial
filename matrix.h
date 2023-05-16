#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

#ifndef CPLUSPLUS_MATRIX_H
#define CPLUSPLUS_MATRIX_H

typedef long long ll;

template <typename T>
struct Rational {
    T p = 0;
    T q = 1;
    Rational(T val) {
        p = val;
        q = 1;
    }
    Rational(T num, T den) {
        p = num;
        q = den;
        fix();
    }
    Rational& fix() {
        if (p == 0 || q == 0) {
            p = 0;
            q = 1;
        }
        else {
            T g = std::__gcd(p, q);
            p /= g;
            q /= g;
            if (q < 0) {
                p = -p;
                q = -q;
            }
        }
        return *this;
    }
    Rational operator-() const {
        return Rational(-p, q);
    }
    friend Rational operator+(const Rational& a, const Rational& b) {
        return Rational(a.p * b.q + a.q * b.p, a.q * b.q);
    }
    friend Rational operator-(const Rational& a, const Rational& b) {
        return Rational(a.p * b.q - a.q * b.p, a.q * b.q);
    }
    friend Rational operator*(const Rational& a, const Rational& b) {
        return Rational(a.p * b.p, a.q * b.q);
    }
    friend Rational operator/(const Rational& a, const Rational& b) {
        return Rational(a.p * b.q, a.q * b.p);
    }
    Rational& operator+=(const Rational& ot) {
        return *this = *this + ot;
    }
    Rational& operator-=(const Rational& ot) {
        return *this = *this - ot;
    }
    Rational& operator*=(const Rational& ot) {
        return *this = *this * ot;
    }
    Rational& operator/=(const Rational& ot) {
        return *this = *this / ot;
    }
    Rational& operator++(int) {
        return *this += 1;
    }
    friend bool operator==(const Rational& a, const Rational& b) {
        return a.p == b.p && a.q == b.q;
    }
    friend bool operator!=(const Rational& a, const Rational& b) {
        return a.p != b.p || a.q != b.q;
    }
    friend bool operator<(const Rational& a, const Rational& b) {
        return a.p * b.q < a.q * b.p;
    }
    friend bool operator>(const Rational& a, const Rational& b) {
        return a.p * b.q > a.q * b.p;
    }
};

template <typename T>
Rational<T> abs(const Rational<T>& r) {
    return {abs(r.p), abs(r.q)};
}

template <typename T = long long>
std::istream& operator>>(std::istream& istr, Rational<T>& r) {
    istr >> r.p;
    r.q = 1;
    return istr;
}

template <typename T = long long>
std::ostream& operator<<(std::ostream& ostr, const Rational<T>& r) {
    if (r.q == 1) {
        return ostr << r.p;
    } else {
        return ostr << "(" << r.p << " / " << r.q << ")";
    }
}

template <typename T>
struct Matrix {
    int row, col;
    std::vector<std::vector<T>> data;
    Matrix(int _row = 0, int _col = 0, bool is_identity = 0) {
        row = _row;
        col = _col;
        data.resize(row);
        for (int i = 0; i < row; i++) {
            data[i].resize(col, 0);
        }
        if (is_identity) {
            for (int i = 0; i < row; i++) {
                data[i][i] = 1;
            }
        }
    }
    std::vector<T>& operator[](size_t i) {
        return data[i];
    }
    const std::vector<T>& operator[](size_t i) const {
        return data[i];
    }
};

template <typename T>
std::istream& operator>>(std::istream& istr, Matrix<T>& a) {
    for (int i = 0; i < a.row; i++) {
        for (int j = 0; j < a.col; j++) {
            istr >> a[i][j];
        }
    }
    return istr;
}

template <typename T>
std::ostream& operator<<(std::ostream& ostr, const Matrix<T>& a) {
    for (int i = 0; i < a.row; i++) {
        if (i != 0) ostr << "\n";
        for (int j = 0; j < a.col; j++) {
            if (j != 0) ostr << " ";
            ostr << a[i][j];
        }
    }
    return ostr;
}

// transforms matrix A to its reduced row echelon form using Gaussian elimination and returns its determinant
template <typename T>
T Gauss(Matrix<T>& a) {
    T det = 1;
    for (int row = 0, col = 0; row < a.row && col < a.col; col++) {
        int sel = row;
        for (int i = row; i < a.row; i++) {
            if (abs(a[i][col]) > abs(a[sel][col])) {
                sel = i;
            }
        }
        if (a[sel][col] == T(0)) {
            continue;
        }
        if (sel != row) {
            det *= T(-1);
        }
        for (int j = col; j < a.col; j++) {
            std::swap(a[sel][j], a[row][j]);
        }
        for (int i = 0; i < a.row; i++) {
            if (i != row) {
                T c = a[i][col] / a[row][col];
                for (int j = col; j < a.col; j++) {
                    a[i][j] -= a[row][j] * c;
                }
            }
        }
        row++;
    }
    for (int i = 0; i < a.row; i++) {
        T val = 0;
        for (int j = 0; j < a.col; j++) {
            if (a[i][j] != T(0)) {
                val = T(1) / a[i][j];
                break;
            }
        }
        if (val == T(0)) {
            det = 0;
        }
        else {
            det /= val;
        }
        for (int j = 0; j < a.col; j++) {
            a[i][j] *= val;
        }
    }
    return det;
}

#define rank my_rank

template <typename T>
int rank(Matrix<T> A) {
    Gauss(A);
    for (int i = 0; i < A.row; i++) {
        bool ok = false;
        for (int j = 0; j < A.col; j++) {
            if (A[i][j] != T(0)) {
                ok = true;
                break;
            }
        }
        if (!ok) {
            return i;
        }
    }
    return A.row;
}


#endif //CPLUSPLUS_MATRIX_H
