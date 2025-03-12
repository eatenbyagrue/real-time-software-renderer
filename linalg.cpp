#include "linalg.hpp"
#include <array>
#include <cassert>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>

Vec4D::Vec4D() { this->set_0(); }

Vec4D::Vec4D(std::array<double, 4> array) : v(array) {}

std::string Vec4D::to_string() {
    return "x: " + std::to_string(this->v.at(0)) + " y: " + std::to_string(this->v.at(1)) +
           " z: " + std::to_string(this->v.at(2)) + " w: " + std::to_string(this->v.at(3));
}

Vec4D Vec4D::operator-() {
    return Vec4D({-this->v.at(0), -this->v.at(1), -this->v.at(2), this->v.at(3)});
}

void Vec4D::operator+=(const Vec4D &rhs) {
    this->v.at(0) += rhs.v.at(0);
    this->v.at(1) += rhs.v.at(1);
    this->v.at(2) += rhs.v.at(2);
    this->v.at(3) += rhs.v.at(3);
}
void Vec4D::operator/=(double rhs) {
    this->v.at(0) /= rhs;
    this->v.at(1) /= rhs;
    this->v.at(2) /= rhs;
    this->v.at(3) /= rhs;
}

Vec4D Vec4D::operator*(const double rhs) {
    return Vec4D({
        this->v.at(0) * rhs,
        this->v.at(1) * rhs,
        this->v.at(2) * rhs,
        this->v.at(3) * rhs,
    });
}
Vec4D Vec4D::operator-(const Vec4D &rhs) {
    return Vec4D({
        this->v.at(0) - rhs.v.at(0),
        this->v.at(1) - rhs.v.at(1),
        this->v.at(2) - rhs.v.at(2),
        this->v.at(3) - rhs.v.at(3),
    });
}
Vec4D Vec4D::operator+(const Vec4D &rhs) {
    return Vec4D({
        this->v.at(0) + rhs.v.at(0),
        this->v.at(1) + rhs.v.at(1),
        this->v.at(2) + rhs.v.at(2),
        this->v.at(3) + rhs.v.at(3),
    });
}

void Vec4D::set_0() {
    this->v.at(0) = 0;
    this->v.at(1) = 0;
    this->v.at(2) = 0;
    this->v.at(3) = 0;
}

double Vec4D::norm() {
    // Only works for vectors!
    assert(this->v.at(3) == 0);
    return sqrt(this->dot(*this));
}

double Vec4D::distance(Vec4D vector) {
    // Assert that both are either vector or point (w coordinate is the same)
    assert(this->v.at(3) == vector.v.at(3));
    Vec4D difference = *this - vector;
    return difference.norm();
}

double Vec4D::dot(Vec4D rhs) {
    // Surprise, just a wrapper for the lengthy std function
    return std::inner_product(this->v.begin(), this->v.end(), rhs.v.begin(), 0.0);
}

Plane::Plane(Vec4D n, double d) : n(n), d(d) {}

void Plane::normalize() {
    double norm = sqrt(this->n.dot(this->n));
    this->n /= norm;
}

double Plane::signed_distance(Vec4D point) { return this->n.dot(point) + this->d; }

Vec4D Plane::intersect(Vec4D A, Vec4D B) {
    Vec4D result;
    double t = (-this->d - this->n.dot(A)) / this->n.dot(B - A);
    result = A + ((B - A) * t);
    return result;
}

Matrix4D::Matrix4D() {
    this->rows = 4;
    this->cols = 4;
    for (int i = 0; i < this->rows * this->cols; ++i) {
        this->m.push_back(0);
    }
};

Matrix4D::Matrix4D(int rows, int cols, std::vector<double> values)
    : rows(rows), cols(cols), m(values) {};

Matrix4D::Matrix4D(int rows, int cols) {
    this->rows = rows;
    this->cols = cols;
    for (int i = 0; i < rows * cols; ++i) {
        this->m.push_back(0);
    }
}

Matrix4D::Matrix4D(int rows, int cols, double value) {
    this->rows = rows;
    this->cols = cols;
    for (int i = 0; i < rows * cols; ++i) {
        this->m.push_back(value);
    }
}

std::string Matrix4D::to_string() {
    std::string str{};
    for (int i = 0; i < this->rows * this->cols; ++i) {
        str += std::to_string(m.at(i)) + " ";
        i % this->cols == this->cols - 1 ? str += ("\n") : str + " ";
    }
    return str;
}

/*This can be optimized by storing the rows and columns redundantly for faster access*/
std::vector<double> Matrix4D::get_row(int r) {
    std::vector<double> result;
    for (int i = r * this->cols; i < (r + 1) * (this->cols); ++i) {
        result.push_back(m.at(i));
    }
    return result;
}

/*This can be optimized by storing the rows and columns redundantly for faster access*/
std::vector<double> Matrix4D::get_col(int c) {
    std::vector<double> result;
    for (int i = c; i < this->rows * this->cols; i += this->cols) {
        result.push_back(m.at(i));
    }
    return result;
}

// Returns transpose
Matrix4D Matrix4D::T() {
    Matrix4D result = Matrix4D();

    // set diagonal
    result.set(0, 0, this->get(0, 0));
    result.set(1, 1, this->get(1, 1));
    result.set(2, 2, this->get(2, 2));
    result.set(3, 3, this->get(3, 3));

    // set bottom left of diagonal
    result.set(0, 1, this->get(1, 0));
    result.set(0, 2, this->get(2, 0));
    result.set(0, 3, this->get(3, 0));
    result.set(1, 2, this->get(2, 1));
    result.set(1, 3, this->get(3, 1));
    result.set(2, 3, this->get(3, 2));

    // set top right of diagonal
    result.set(1, 0, this->get(0, 1));
    result.set(2, 0, this->get(0, 2));
    result.set(3, 0, this->get(0, 3));
    result.set(2, 1, this->get(1, 2));
    result.set(3, 1, this->get(1, 3));
    result.set(3, 2, this->get(2, 3));

    return result;
}

double Matrix4D::get(int row, int col) { return m.at(row * cols + col); }

void Matrix4D::set(int row, int col, double value) { m.at(row * cols + col) = value; }

void Matrix4D::set_3D(std::array<double, 9> values) {
    for (int i = 0; i < values.size(); ++i) {
        m.at(i) = values.at(i);
    }
}

// Overload the << operator
std::ostream &operator<<(std::ostream &os, Matrix4D &obj) {
    os << obj.to_string(); // Print the internal string
    return os;
}

Matrix4D Matrix4D::operator*(Matrix4D &rhs) {
    assert(this->cols == rhs.rows);
    Matrix4D result(this->rows, rhs.cols);
    for (int r = 0; r < this->rows; ++r) {
        for (int c = 0; c < rhs.cols; ++c) {
            std::vector get_row = this->get_row(r);
            double dotpr =
                std::inner_product(get_row.begin(), get_row.end(), rhs.get_col(c).begin(), 0.0);
            result.set(r, c, dotpr);
        }
    }
    return result;
}

Vec4D Matrix4D::times(Vec4D &hovec) {
    assert(this->cols == hovec.v.size());
    assert(this->rows <= hovec.v.size());
    Vec4D result;
    for (int r = 0; r < this->rows; ++r) {
        double sum = 0.0;
        for (int c = 0; c < this->cols; ++c) {
            sum += this->get(r, c) * hovec.v.at(c);
        }
        result.v.at(r) = sum;
        // For some reason the dot product doesnt give the right result here. what's going wrong?
        /*std::vector get_row = this->get_row(r);*/
        /*result.v.at(r) = std::inner_product(hovec.v.begin(), hovec.v.end(), get_row.begin(),
         * 0.0);*/
    }
    return result;
}

Matrix4D Matrix4D::create_identity_matrix(int rows) {
    Matrix4D id{rows, rows};
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            i == j ? id.set(i, j, 1.0) : id.set(i, j, 0.0);
        }
    }
    return id;
}

/* Rotation Matrix4D in 3D around the axis, by alpha degrees.
 * 0 for x axis, 1 for y axis, 2 for z axis */
Matrix4D Matrix4D::create_rotation_matrix(double alpha, int axis) {
    std::vector<double> values;
    switch (axis) {
    case 0:
        values = {1, 0,          0,          0, 0, cos(alpha), -sin(alpha), 0,
                  0, sin(alpha), cos(alpha), 0, 0, 0,          0,           1};
        break;
    case 1:
        values = {cos(alpha),  0, sin(alpha), 0, 0, 1, 0, 0,
                  -sin(alpha), 0, cos(alpha), 0, 0, 0, 0, 1};
        break;
    case 2:
        values = {cos(alpha), -sin(alpha), 0, 0, sin(alpha), cos(alpha), 0, 0,
                  0,          0,           1, 0, 0,          0,          0, 1};
        break;
    }
    Matrix4D rot{4, 4, values};
    return rot;
}

Matrix4D Matrix4D::create_scale_matrix(double scale) {
    Matrix4D sc{4, 4, 0.0};
    sc.set(0, 0, scale);
    sc.set(1, 1, scale);
    sc.set(2, 2, scale);
    sc.set(3, 3, 1.0);
    return sc;
}

Matrix4D Matrix4D::create_translation_matrix(Vec4D translation) {
    return Matrix4D::create_translation_matrix(
        {translation.v.at(0), translation.v.at(1), translation.v.at(2)});
};
Matrix4D Matrix4D::create_translation_matrix(std::array<double, 3> translation) {
    Matrix4D tr{4, 4, 0.0};
    tr.set(0, 0, 1.0);
    tr.set(1, 1, 1.0);
    tr.set(2, 2, 1.0);
    tr.set(3, 3, 1.0);
    tr.set(0, 3, translation.at(0));
    tr.set(1, 3, translation.at(1));
    tr.set(2, 3, translation.at(2));
    return tr;
}
