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

//template<uword N, typename T>
//std::array<T, N> vec_to_array()
