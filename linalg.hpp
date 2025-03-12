#ifndef LINALG_H
#define LINALG_H

#include <array>
#include <string>
#include <vector>

class Vec3D {
    std::array<double, 3> v;
};

/*A 4-dimensional vector in homogenous coordinates*/
class Vec4D {
  public:
    std::array<double, 4> v;
    Vec4D();
    Vec4D(std::array<double, 4> array);
    std::string to_string();
    Vec4D operator-();
    void operator+=(const Vec4D &rhs);
    Vec4D operator-(const Vec4D &rhs);
    Vec4D operator+(const Vec4D &rhs);
    void operator/=(double rhs);
    Vec4D operator*(const double rhs);
    /*Sets all values of v to 0*/
    void set_0();
    /*Calculates the length of the vector.*/
    /*Only works for vectors, i.e. v.at(3) == 0. */
    double norm();
    double distance(Vec4D vector);
    /*Calculates the dot product between this vector and another*/
    /*Careful use required if 4th coordinate is not equal to 0.0, since it is used in the
     * calculation, too!*/
    double dot(Vec4D rhs);
};

class Plane {
  private:
    Vec4D n;
    double d;

  public:
    Plane(Vec4D n, double d);
    // normalizes the plane normal to unit length
    void normalize();
    // Calculates the distance of the point from the plane
    double signed_distance(Vec4D point);
    // Calculates the intersection of the plane and a line given by two points
    Vec4D intersect(Vec4D A, Vec4D B);
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
