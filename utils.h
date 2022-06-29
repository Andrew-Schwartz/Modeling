#include <armadillo>

using namespace arma;

/// for parsing utilites
struct str {
  const std::string &string;

  explicit str(const std::string &string) : string(string) {}

  bool starts_with(const std::string &needle) {
    return string.find(needle) == 0;
  }

  bool contains(const std::string &needle) {
    return string.find(needle) != std::string::npos;
  }

//  const std::string &operator*() {
//    return string;
//  }

  const std::string *operator->() {
    return &string;
  }
};

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

// todo make this work for any T - idk why it didn't
//template<typename T>
std::optional<double> map_opt(std::optional<double> opt, const std::function<double(double)> &f) {
  return opt
         ? std::make_optional(f(*opt))
         : std::nullopt;
}

//template<uword N, typename T>
//std::array<T, N> vec_to_array()