#include "SDL2/SDL.h"
/*#include "SDL_events.h"*/
/*#include "SDL_render.h"*/
#include "SDL_scancode.h"
#include "linalg.hpp"
#include "profiler.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h>
#include <ostream>
#include <random>
#include <stdio.h>
#include <vector>

/*#define WINDOW_WIDTH 640*/
/*#define WINDOW_HEIGHT 480*/
#define WINDOW_WIDTH 1920.0
#define WINDOW_HEIGHT 1080.0
/*#define VIEWPORT_WIDTH 16.0)*/
/*#define VIEWPORT_HEIGHT 9*/
#define VIEWPORT_WIDTH 8.0
#define VIEWPORT_HEIGHT 4.5
#define DISTANCE_D 4.5
#define PI 3.14159265

struct Color {
    int r, g, b, a;
};

enum class Visibility { FULL, PARTIAL, NONE };

Color RED{255, 0, 0, 0};
Color GREEN{0, 255, 0, 0};
Color BLUE{0, 0, 255, 0};
Color YELLOW{255, 255, 0, 0};
Color PURPLE{128, 0, 128, 0};
Color CYAN{0, 255, 255, 0};
Color WHITE{255, 255, 255, 0};

enum class ModelType {
    CUBE,
    SPHERE,
    CUBOID,
    CARTESIAN,
};

struct Point {
    // x,y are 2-D coordinates
    // h is the shading, 0.0 for black and 1.0 for normal color
    int x, y;
    double h;

    void swap(Point &p) {
        Point _{this->x, this->y, this->h};
        this->x = p.x;
        this->y = p.y;
        this->h = p.h;
        p.x = _.x;
        p.y = _.y;
        p.h = _.h;
    }

    Point flip() { return Point{this->y, this->x, this->h}; }

    std::string to_string() {
        return "x: " + std::to_string(this->x) + " y: " + std::to_string(this->y);
    }
};

struct Triangle {
    int A, B, C;
    Color color;
    Triangle(int A, int B, int C, Color color) : A(A), B(B), C(C), color(color) {};
};

struct Model {
    ModelType model_type;
    std::vector<Vec4D> vertices;
    std::vector<Triangle> triangles;
};

struct Transform {
    double scale{1};
    double rotation_y{0.0};
    std::array<double, 3> translation;
    /*Transform(double s, double r, std::array<double, 3> l)*/
    /*: scale(s), rotation_y(r), translation(l) {};*/
};

class Instance {
  public:
    std::vector<Vec4D> vertices;
    std::vector<Triangle> triangles;
    Vec4D bounding_sphere_center;
    double bounding_sphere_radius;
    Model model;
    Transform transform;
    Instance(Model model, Transform transform) : model(model), transform(transform) {};

    // Averages all vertices
    // This is terribly inefficient to do each frame
    Vec4D calc_center() {
        Vec4D result;
        for (Vec4D &vertex : this->vertices) {
            result += vertex;
        }
        result /= this->vertices.size();
        return result;
    }

    // Find the distance to the vertex furthest away from the center
    double calc_radius(Vec4D center) {
        double result{0};
        for (Vec4D &vertex : this->vertices) {
            double dist = vertex.distance(center);
            if (dist > result)
                result = dist;
        }
        return result;
    }
};

struct Camera {
    Vec4D position;
    Vec4D orientation;
    Vec4D original_orientation;
    // This is just a crutch until I figure out full 3d movement of the camera
    double rotation_y;
    std::vector<Plane> planes;
};

struct Scene {
    std::vector<Model> models;
    std::vector<Instance> instances;
    Camera camera;
};

std::vector<double> interpolate(int i0, double d0, int i1, double d1) {
    if (i0 == i1) {
        return {d0};
    }
    std::vector<double> values{};
    /*values.reserve(1 + i1 - i0);*/

    double a{((double)d1 - d0) / ((double)i1 - i0)};
    double d = d0;
    for (int i = i0; i <= i1; ++i) {
        values.push_back(d);
        d = d + a;
    }
    return values;
}

/*Draw point with x and y coordinates*/
void draw_point(SDL_Renderer *renderer, int x, int y) {
    if (-WINDOW_WIDTH / 2 > x || x >= WINDOW_WIDTH / 2 || -WINDOW_HEIGHT / 2 > y ||
        y >= WINDOW_HEIGHT / 2) {
        return;
    }
    int draw_x = WINDOW_WIDTH / 2 + x;
    int draw_y = WINDOW_HEIGHT / 2 - y;
    SDL_RenderDrawPoint(renderer, draw_x, draw_y);
}

/*void draw_point(SDL_Renderer *renderer, int x, int y, Color &color) {*/
/*    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);*/
/*    SDL_RenderDrawPoint(renderer, x, y);*/
/*}*/

void draw_line(SDL_Renderer *renderer, Point p0, Point p1) {
    int dx{p1.x - p0.x};
    int dy{p1.y - p0.y};

    if (abs(dx) > abs(dy)) {
        // Line is somewhat horizontal-ish
        if (p0.x > p1.x) {
            p0.swap(p1);
        }
        std::vector ys = interpolate(p0.x, p0.y, p1.x, p1.y);
        for (int x = p0.x; x <= p1.x; ++x) {
            draw_point(renderer, x, round(ys.at(x - p0.x)));
        }
    } else {
        // Line is somewhat vertical-ish
        if (p0.y > p1.y) {
            p0.swap(p1);
        }
        std::vector xs = interpolate(p0.y, p0.x, p1.y, p1.x);
        for (int y = p0.y; y <= p1.y; ++y) {
            draw_point(renderer, round(xs.at(y - p0.y)), y);
        }
    }
}

std::vector<Point> gen_random_points(int count) {
    // random
    std::vector<Point> random_points;
    /*std::random_device random_device;*/
    std::mt19937 generator(31337);
    /*std::mt19937 generator(random_device());*/
    for (int i = 0; i < count; ++i) {
        std::uniform_int_distribution<int> dist_x(-WINDOW_WIDTH / 2, WINDOW_WIDTH / 2 - 1);
        std::uniform_int_distribution<int> dist_y(-WINDOW_HEIGHT / 2, WINDOW_HEIGHT / 2 - 1);
        std::uniform_real_distribution<double> dist_h(0.0, 1.0);
        Point point{dist_x(generator), dist_y(generator), dist_h(generator)};
        random_points.push_back(point);
        /*std::cout << point.x << " " << point.y << " " << point.h << std::endl;*/
    }
    return random_points;
}

void draw_wireframe_triangle(SDL_Renderer *renderer, Point p0, Point p1, Point p2) {
    draw_line(renderer, p0, p1);
    draw_line(renderer, p1, p2);
    draw_line(renderer, p2, p0);
}

void draw_wireframe_triangle(SDL_Renderer *renderer, Point p0, Point p1, Point p2, Color color) {
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
    draw_line(renderer, p0, p1);
    draw_line(renderer, p1, p2);
    draw_line(renderer, p2, p0);
}

std::vector<double> concatenate(std::vector<double> &vec0, std::vector<double> &vec1) {
    std::vector<double> concat;
    concat.reserve(vec0.size() + vec1.size());
    // concatenate vec0 and vec1 into concat
    concat.insert(concat.end(), vec0.begin(), vec0.end());
    concat.insert(concat.end(), vec1.begin(), vec1.end());
    return concat;
}

void draw_filled_triangle(SDL_Renderer *renderer, Point p0, Point p1, Point p2, Color color) {
    if (p1.y < p0.y)
        p1.swap(p0);
    if (p2.y < p0.y)
        p2.swap(p0);
    if (p2.y < p1.y)
        p2.swap(p1);

    std::vector<double> x01 = interpolate(p0.y, p0.x, p1.y, p1.x);
    std::vector<double> x12 = interpolate(p1.y, p1.x, p2.y, p2.x);
    std::vector<double> x02 = interpolate(p0.y, p0.x, p2.y, p2.x);

    // remove last element of x01, since it is duplucated in x12
    x01.pop_back();
    std::vector<double> x012 = concatenate(x01, x12);

    int m = std::floor(x012.size() / 2);
    if (x02.at(m) < x012.at(m)) {
        for (int y = p0.y; y <= p2.y; ++y) {
            for (int x = round(x02[y - p0.y]); x <= x012[y - p0.y]; ++x) {
                SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
                draw_point(renderer, x, y);
            }
        }
    } else {
        for (int y = p0.y; y <= p2.y; ++y) {
            for (int x = round(x012[y - p0.y]); x <= x02[y - p0.y]; ++x) {
                SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
                draw_point(renderer, x, y);
            }
        }
    }
}

/* There are some bugs in here! */
void draw_shaded_triangle(SDL_Renderer *renderer, Point &p0, Point &p1, Point &p2, Color &color) {
    if (p1.y < p0.y)
        p1.swap(p0);
    if (p2.y < p0.y)
        p2.swap(p0);
    if (p2.y < p1.y)
        p2.swap(p1);

    std::vector<double> x01 = interpolate(p0.y, p0.x, p1.y, p1.x);
    std::vector<double> h01 = interpolate(p0.y, p0.h, p1.y, p1.h);
    std::vector<double> x12 = interpolate(p1.y, p1.x, p2.y, p2.x);
    std::vector<double> h12 = interpolate(p1.y, p1.h, p2.y, p2.h);
    std::vector<double> x02 = interpolate(p0.y, p0.x, p2.y, p2.x);
    std::vector<double> h02 = interpolate(p0.y, p0.h, p2.y, p2.h);

    x01.pop_back();
    std::vector<double> x012 = concatenate(x01, x12);
    h01.pop_back();
    std::vector<double> h012 = concatenate(h01, h12);

    // inefficient, change eventually
    std::vector<double> x_left;
    std::vector<double> h_left;
    std::vector<double> x_right;
    std::vector<double> h_right;

    int m = std::floor(x012.size() / 2);
    if (x02.at(m) < x012.at(m)) {
        x_left = x02;
        h_left = h02;
        x_right = x012;
        h_right = h012;
    } else {
        x_left = x012;
        h_left = h012;
        x_right = x02;
        h_right = h02;
    }
    for (int y = p0.y; y <= p2.y; ++y) {
        double x_l = x_left[y - p0.y];
        double x_r = x_right[y - p0.y];
        std::vector<double> h_segment = interpolate(x_l, h_left[y - p0.y], x_r, h_right[y - p0.y]);
        for (int x = x_l; x <= x_r; ++x) {
            double h = h_segment[x - x_l];
            SDL_SetRenderDrawColor(renderer, h * color.r, h * color.g, h * color.b, 0);
            draw_point(renderer, x, y);
        }
    }
}

Point viewport_to_canvas(double x, double y) {
    return Point{static_cast<int>(x * WINDOW_WIDTH / VIEWPORT_WIDTH),
                 static_cast<int>(y * WINDOW_HEIGHT / VIEWPORT_HEIGHT), 1.0};
}

Point project_vertex(Vec4D &h) {
    return viewport_to_canvas(h.v.at(0) * DISTANCE_D / h.v.at(2),
                              h.v.at(1) * DISTANCE_D / h.v.at(2));
}

Model create_cube() {
    Model model;
    model.model_type = ModelType::CUBE;
    std::array<double, 4> derp{1, 1, 1, 1};
    model.vertices.push_back(Vec4D({1.0, 1.0, 1.0, 1.0}));
    model.vertices.push_back(Vec4D({-1, 1, 1, 1}));
    model.vertices.push_back(Vec4D({-1, -1, 1, 1}));
    model.vertices.push_back(Vec4D({1, -1, 1, 1}));
    model.vertices.push_back(Vec4D({1, 1, -1, 1}));
    model.vertices.push_back(Vec4D({-1, 1, -1, 1}));
    model.vertices.push_back(Vec4D({-1, -1, -1, 1}));
    model.vertices.push_back(Vec4D({1, -1, -1, 1}));

    model.triangles.emplace_back(0, 1, 2, RED);
    model.triangles.emplace_back(0, 2, 3, RED);
    model.triangles.emplace_back(4, 0, 3, GREEN);
    model.triangles.emplace_back(4, 3, 7, GREEN);
    model.triangles.emplace_back(5, 4, 7, BLUE);
    model.triangles.emplace_back(5, 7, 6, BLUE);
    model.triangles.emplace_back(1, 5, 6, YELLOW);
    model.triangles.emplace_back(1, 6, 2, YELLOW);
    model.triangles.emplace_back(4, 5, 1, PURPLE);
    model.triangles.emplace_back(4, 1, 0, PURPLE);
    model.triangles.emplace_back(2, 6, 7, CYAN);
    model.triangles.emplace_back(2, 7, 3, CYAN);

    return model;
};

Scene create_scene() {
    Scene scene;
    Camera camera;
    camera.position = Vec4D({0, 0, 0, 1});
    camera.orientation = Vec4D({0, 0, 1, 0});
    camera.original_orientation = Vec4D({0, 0, 1, 0});
    camera.rotation_y = 0.0;

    // Define camera frustrum
    std::vector<Plane> planes;
    Plane p_near{Vec4D({0, 0, 1, 0}), -DISTANCE_D / 2};
    p_near.normalize();
    planes.push_back(p_near);
    Plane p_left{Vec4D({DISTANCE_D, 0, VIEWPORT_WIDTH / 2.0, 0}), 0};
    p_left.normalize();
    planes.push_back(p_left);
    Plane p_right{Vec4D({-DISTANCE_D, 0, VIEWPORT_WIDTH / 2.0, 0}), 0};
    p_right.normalize();
    planes.push_back(p_right);
    Plane p_bottom{Vec4D({0, DISTANCE_D, VIEWPORT_HEIGHT / 2.0, 0}), 0};
    p_bottom.normalize();
    planes.push_back(p_bottom);
    Plane p_top{Vec4D({0, -DISTANCE_D, VIEWPORT_HEIGHT / 2.0, 0}), 0};
    p_top.normalize();

    camera.planes = planes;
    scene.camera = camera;

    Model model = create_cube();

    scene.instances.emplace_back(model, Transform{1.0, 0.3, {-2, 0, 5.0}});
    scene.instances.emplace_back(model, Transform{0.5, PI, {-2, 1, 8.0}});
    scene.instances.emplace_back(model, Transform{1.0, 0.3, {5, -1, 10.0}});

    return scene;
}

// Change the 2nd parameter to camera.orientation once I figure out camera movement
Matrix4D make_camera_matrix(Vec4D position, double rotation_y) {
    Matrix4D m_rotation = Matrix4D::create_rotation_matrix(rotation_y, 1);
    Matrix4D m_translation = Matrix4D::create_translation_matrix(-position);
    return m_rotation * m_translation;
}

Matrix4D make_instance_transform_matrix(Transform transform) {
    Matrix4D m_scale = Matrix4D::create_scale_matrix(transform.scale);
    Matrix4D m_rotation = Matrix4D::create_rotation_matrix(transform.rotation_y, 1);
    Matrix4D m_translation = Matrix4D::create_translation_matrix(transform.translation);
    return m_translation * m_rotation * m_scale;
}

Instance transform_instance(Instance &instance, Matrix4D transform) {
    Instance transformed_instance = instance;
    for (Vec4D &vertex : transformed_instance.model.vertices) {
        Vec4D transformed_vertex = transform.times(vertex);
        transformed_instance.vertices.push_back(transformed_vertex);
    }
    transformed_instance.triangles = instance.model.triangles;
    return transformed_instance;
}

Instance clip_instance_against_plane(Instance &instance, Plane &plane) {
    Instance clipped_instance = instance;
    std::vector<Triangle> clipped_triangles;

    for (Triangle &triangle : instance.triangles) {
        // save the distances and triangle point indices in pairs
        std::array<std::pair<double, int>, 3> distances_vertices;
        distances_vertices.at(0) = {plane.signed_distance(clipped_instance.vertices.at(triangle.A)),
                                    triangle.A};
        distances_vertices.at(1) = {plane.signed_distance(clipped_instance.vertices.at(triangle.B)),
                                    triangle.B};
        distances_vertices.at(2) = {plane.signed_distance(clipped_instance.vertices.at(triangle.C)),
                                    triangle.C};

        // Sort in ascending order by distance
        std::sort(distances_vertices.begin(), distances_vertices.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
        // Count positive distances
        int count_positive = std::count_if(distances_vertices.begin(), distances_vertices.end(),
                                           [](const auto &d) { return d.first >= 0; });

        switch (count_positive) {
        case 3: {
            // If all three distances are positive, keep the triangle as is
            clipped_triangles.push_back(triangle);
            break;
        }
        case 0: {
            // If all three distances are negative, discard the triangle
            break;
        }
        case 1: {
            // If only one distance is positive, make new triangle
            // The positive distance is the last element
            int A_index = distances_vertices.at(2).second;
            int B_index = distances_vertices.at(0).second;
            int C_index = distances_vertices.at(1).second;
            Vec4D A = clipped_instance.vertices.at(A_index);
            Vec4D B = clipped_instance.vertices.at(B_index);
            Vec4D C = clipped_instance.vertices.at(C_index);
            Vec4D B_prime = plane.intersect(A, B);
            Vec4D C_prime = plane.intersect(A, C);
            clipped_instance.vertices.push_back(B_prime);
            int B_prime_index = clipped_instance.vertices.size() - 1;
            clipped_instance.vertices.push_back(C_prime);
            int C_prime_index = clipped_instance.vertices.size() - 1;
            /*Triangle tria =*/
            /*    Triangle(distances_vertices.at(2).second, B_index, C_index, triangle.color);*/
            Triangle tria = Triangle(A_index, B_index, C_index, triangle.color);
            clipped_triangles.push_back(tria);
            break;
        }
        case 2: {
            // If two distances are positive, make two new triangles
            // The negative distance is the first element
            int C_index = distances_vertices.at(0).second;
            int A_index = distances_vertices.at(1).second;
            int B_index = distances_vertices.at(2).second;
            Vec4D C = clipped_instance.vertices.at(C_index);
            Vec4D A = clipped_instance.vertices.at(A_index);
            Vec4D B = clipped_instance.vertices.at(B_index);
            Vec4D A_prime = plane.intersect(A, C);
            Vec4D B_prime = plane.intersect(B, C);
            clipped_instance.vertices.push_back(A_prime);
            int A_prime_index = clipped_instance.vertices.size() - 1;
            clipped_instance.vertices.push_back(B_prime);
            int B_prime_index = clipped_instance.vertices.size() - 1;
            Triangle tria = Triangle(A_index, B_index, A_prime_index, triangle.color);
            Triangle tria2 = Triangle(A_prime_index, B_index, B_prime_index, triangle.color);
            clipped_triangles.push_back(tria);
            clipped_triangles.push_back(tria2);
            break;
        }
        }
    }
    clipped_instance.triangles = clipped_triangles;
    return clipped_instance;
}

Instance clip_instance(Instance &instance, std::vector<Plane> planes) {
    Instance clipped_instance = instance;
    for (Plane &plane : planes) {
        clipped_instance = clip_instance_against_plane(clipped_instance, plane);
    }
    return clipped_instance;
};

void draw_instance(SDL_Renderer *renderer, Instance &instance) {
    std::vector<Point> projected_vertices;

    for (Vec4D &vertex : instance.vertices) {
        projected_vertices.push_back(project_vertex(vertex));
    }
    for (Triangle &triangle : instance.triangles) {
        draw_wireframe_triangle(renderer, projected_vertices.at(triangle.A),
                                projected_vertices.at(triangle.B),
                                projected_vertices.at(triangle.C), triangle.color);
        /*draw_filled_triangle(renderer, projected_vertices.at(triangle.A),*/
        /*                     projected_vertices.at(triangle.B),
         * projected_vertices.at(triangle.C),*/
        /*                     triangle.color);*/
    }
}

Visibility is_instance_visible(Instance &instance, std::vector<Plane> planes) {
    bool partiality{false};
    for (Plane &plane : planes) {
        // Inefficient to do every frame, calculate once on creation and then
        // transform center along with rest
        Vec4D instance_center = instance.calc_center();
        double r = instance.calc_radius(instance_center);
        double d = plane.signed_distance(instance_center);
        if (d < -r)
            // Instance is completely out of view
            return Visibility::NONE;
        else if (d < r) {
            // Instance is partially out of view
            partiality = true;
        }
    }
    // Instance is at least partially in view
    if (partiality)
        return Visibility::PARTIAL;
    return Visibility::FULL;
}

void render_scene(SDL_Renderer *renderer, Scene &scene, Profiler profiler) {
    // clip_instance and draw_instance operate on instance in place

    profiler.start(Profile::TRANSFORMING);
    // transform scene
    Matrix4D m_camera = make_camera_matrix(scene.camera.position, scene.camera.rotation_y);

    std::vector<Instance> transformed_instances;
    for (Instance &instance : scene.instances) {
        Matrix4D m_transform = make_instance_transform_matrix(instance.transform);
        Matrix4D m_f = m_camera * m_transform;
        transformed_instances.push_back(transform_instance(instance, m_f));
    }
    profiler.stop();

    profiler.start(Profile::CLIPPING);
    // clip scene
    std::vector<Instance> clipped_instances;
    for (Instance &instance : transformed_instances) {
        // If instance is fully visible, take it as is.
        // If instance is only partially visible, determine on triangle basis.
        // If instance is completely out of view, discard.
        switch (is_instance_visible(instance, scene.camera.planes)) {
        case (Visibility::FULL):
            clipped_instances.push_back(instance);
            break;
        case (Visibility::PARTIAL):
            clipped_instances.push_back(clip_instance(instance, scene.camera.planes));
            break;
        case (Visibility::NONE):
            break;
        default:
            break;
        }
    }
    profiler.stop();

    profiler.start(Profile::DRAWING);
    // draw scene
    for (Instance &instance : clipped_instances) {
        draw_instance(renderer, instance);
    }
    profiler.stop();

    profiler.draw_profiler_card();
}

int main(int argc, char *argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "Could not init SDL: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *screen = SDL_CreateWindow("My application", SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED, WINDOW_WIDTH, WINDOW_HEIGHT, 0);
    if (!screen) {
        fprintf(stderr, "Could not create window\n");
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(screen, -1, SDL_RENDERER_SOFTWARE);
    if (!renderer) {
        fprintf(stderr, "Could not create renderer\n");
        return 1;
    }

    SDL_Event evt;

    Scene scene = create_scene();
    Profiler profiler{};
    profiler.register_profile(Profile::TRANSFORMING, "Time taken for transforming (ms): ");
    profiler.register_profile(Profile::CLIPPING, "Time taken for clipping (ms): ");
    profiler.register_profile(Profile::DRAWING, "Time taken for drawing (ms): ");

    bool run = true;

    double t{0.0};
    double x_rot_speed = 0.0;
    double y_rot_speed = 0.0;
    double z_offset = 0.0;

    Vec4D camera_speed;
    camera_speed.set_0();
    // game loop
    while (run) {

        // process OS events
        while (SDL_PollEvent(&evt) != 0) {
            switch (evt.type) {
            case SDL_QUIT:
                run = false;
                break;
            case SDL_KEYDOWN:
                if (evt.key.keysym.scancode == SDL_SCANCODE_LEFT) {
                    y_rot_speed = 0.005;
                } else if (evt.key.keysym.scancode == SDL_SCANCODE_RIGHT) {
                    y_rot_speed = -0.005;
                }
                if (evt.key.keysym.scancode == SDL_SCANCODE_UP) {
                    camera_speed = scene.camera.orientation * 0.01;
                } else if (evt.key.keysym.scancode == SDL_SCANCODE_DOWN) {
                    camera_speed = -scene.camera.orientation * 0.01;
                }

                /*printf("sym %i scancode %i\n", evt.key.keysym.sym,
                 * evt.key.keysym.scancode);*/
                /*std::cout << "Scancode: " << evt.key.keysym.scancode << " ("*/
                /*<<
                 * SDL_GetScancodeName(static_cast<SDL_Scancode>(evt.key.keysym.scancode))*/
                /*<< ")\n";*/
                /*std::cout << "Scancode " << event.key.keysym.scancode <<
                 * "\n";*/
                /*std::cout << "Sym " << SDL_GetKeyName(event.key.keysym.sym) <<
                 * "\n";*/
                /*}*/
                break;
            case SDL_KEYUP:
                if (evt.key.keysym.scancode == SDL_SCANCODE_LEFT) {
                    y_rot_speed = 0.0;
                } else if (evt.key.keysym.scancode == SDL_SCANCODE_RIGHT) {
                    y_rot_speed = 0.0;
                }
                if (evt.key.keysym.scancode == SDL_SCANCODE_UP) {
                    camera_speed.set_0();
                } else if (evt.key.keysym.scancode == SDL_SCANCODE_DOWN) {
                    camera_speed.set_0();
                }
                break;
            }
        }

        scene.camera.rotation_y += y_rot_speed;
        // This is a crutch and should be changed once I figure out camera
        // movement!
        scene.camera.orientation = (Matrix4D::create_rotation_matrix(scene.camera.rotation_y, 1))
                                       .T()
                                       .times(scene.camera.original_orientation);
        scene.camera.position += camera_speed;

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
        SDL_RenderClear(renderer);

        /*scene.instances.at(0).transform.rotation_y += y_rot_speed;*/
        /*scene.instances.at(1).transform.scale = 1.5 + sin(t);*/
        /*scene.instances.at(1).update_transform_matrix();*/
        render_scene(renderer, scene, profiler);
        SDL_RenderPresent(renderer);
        t += 0.011;
    }
    // clean up
    SDL_DestroyWindow(screen);
    SDL_Quit();
    return 0;
}
