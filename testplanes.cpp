#include "linalg.hpp"
#include <cmath>
#include <iostream>
#include <ostream>

#define VIEWPORT_WIDTH 8.0
#define VIEWPORT_HEIGHT 4.5
#define DISTANCE_D 4.0
#define PI 3.14159265

int main() {

    std::vector<Plane> planes;
    Plane p_near{{0, 0, 1}, -DISTANCE_D};
    p_near.normalize();
    planes.push_back(p_near);
    Plane p_left{{DISTANCE_D, 0, VIEWPORT_WIDTH / 2.0}, 0};
    p_left.normalize();
    planes.push_back(p_left);
    Plane p_right{{-DISTANCE_D, 0, VIEWPORT_WIDTH / 2.0}, 0};
    p_right.normalize();
    planes.push_back(p_right);
    Plane p_bottom{{0, DISTANCE_D, VIEWPORT_HEIGHT / 2.0}, 0};
    p_bottom.normalize();
    planes.push_back(p_bottom);
    Plane p_top{{0, -DISTANCE_D, VIEWPORT_HEIGHT / 2.0}, 0};
    p_top.normalize();
    planes.push_back(p_top);
    for (Plane &plane : planes) {
        Vec4D point = Vec4D({DISTANCE_D, 0, VIEWPORT_WIDTH / 2.0, 0});
        std::cout << plane.signed_distance(point) << std::endl;
    }
}
