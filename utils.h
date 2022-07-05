#include <armadillo>

using namespace arma;

template<uword N, typename T>
typename Col<T>::template fixed<N> array_to_vec(const std::array<T, N> &array) {
  typename Col<T>::template fixed<N> col;
  for (int i = 0; i < N; ++i) {
    col(i) = array[i];
  }
  return col;
}

template<uword N, typename T>
bool eq(typename Col<T>::template fixed<N> a, typename Col<T>::template fixed<N> b) {
  for (int i = 0; i < N; ++i)
    if (a(i) != b(i))
      return false;
  return true;
}

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
vec3 rotate_vector(const vec3 &axis, double theta, const vec3 &input) {
  return axis * (dot(axis, input)) + cos(theta) * cross(cross(axis, input), axis) + sin(theta) * cross(axis, input);
}

template<typename ...T>
std::string join_string(const T&... args) {
  std::ostringstream stream;
  (stream << ... << args);
  return stream.str();
}