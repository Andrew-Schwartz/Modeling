#include <armadillo>

using namespace arma;

vec3 rotate_vector(const vec3 &axis, double theta, const vec3 &input) {
    return axis * (dot(axis, input)) + cos(theta) * cross(cross(axis, input), axis) + sin(theta) * cross(axis, input);
}