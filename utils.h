#include <armadillo>
#include <thread>

using namespace arma;

template<uword N, typename T, typename R>
typename Col<R>::template fixed<N> conv_fixed(typename Col<T>::template fixed<N> input) {
    typename Col<R>::template fixed<N> ret;
    for (uword i = 0; i < N; ++i)
        ret(i) = static_cast<R>(input(i));
    return ret;
}

/// Rotates `input` by `theta` around `axis`.
///
/// Rotation math from end of https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
vec3 rotate_vector(const vec3 &axis, double theta, const vec3 &input);

template<typename... Ts>
std::string join_string(const Ts &... args) {
    std::ostringstream stream;
    (stream << ... << args);
    return stream.str();
}

namespace timer {

using namespace std::chrono;
using namespace std::chrono_literals;

template<typename F>
double time(F f) {
    const auto start = high_resolution_clock::now();
    f();
    const auto elapsed = high_resolution_clock::now() - start;
    const double seconds = (double) duration_cast<nanoseconds>(elapsed).count() / 1e9;

    std::cout << "elapsed = " << seconds << std::endl;

    return seconds;
}

}