#ifndef LINALG_H
#define LINALG_H

#include <array>
#include <string>
#include <vector>

// Just use a Vec4D later on
class Vertex {
  public:
    double x, y, z;
    Vertex(double x, double y, double z);
    std::string to_string();
    /* Returns a new vertex. */
    Vertex add(double rhs);
    /* Returns a new vertex */
    Vertex add(std::array<double, 3> rhs);
};

/*A 4-dimensional vector in homogenous coordinates*/
class Vec4D {
  public:
    std::array<double, 4> v;
    Vec4D();
    Vec4D(std::array<double, 4> array);
    Vec4D(Vertex &vertex);
    std::string to_string();
    Vec4D operator-();
    void operator+=(const Vec4D &rhs);
    Vec4D operator*(const double rhs);
    void set_0();
};

class Matrix4D {
    std::vector<double> m;
    int rows, cols;

  public:
    Matrix4D();
    Matrix4D(int rows, int cols, std::vector<double> values);
    Matrix4D(int rows, int cols);
    Matrix4D(int rows, int cols, double value);
    std::string to_string();
    std::vector<double> get_row(int r);
    std::vector<double> get_col(int c);

    Matrix4D T();
    double get(int row, int col);
    void set(int row, int col, double value);
    void set_3D(std::array<double, 9> values);

    friend std::ostream &operator<<(std::ostream &os, Matrix4D &obj);
    Matrix4D operator*(Matrix4D &rhs);
    Vec4D times(Vec4D &rhs);
    static Matrix4D create_identity_matrix(int rows);
    static Matrix4D create_scale_matrix(double scale);
    static Matrix4D create_rotation_matrix(double alpha, int axis);
    static Matrix4D create_translation_matrix(Vec4D translation);
    static Matrix4D create_translation_matrix(std::array<double, 3> translation);
};

#endif
