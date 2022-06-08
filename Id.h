#ifndef MODELING_ID_H
#define MODELING_ID_H

#include <armadillo>

using namespace arma;

constexpr uword SizeX = 41; // 20;
constexpr uword SizeY = 41; // 20;
constexpr uword SizeZ = 41; // 40;
constexpr uword Sizes[3]{SizeX, SizeY, SizeZ};

struct Id {
  uword id;

  explicit Id(uword id);

  Id(uword x, uword y, uword z);

  static Id Invalid;

  [[nodiscard]] bool is_valid() const;

  Id operator+(uword other) const;
  Id operator+(Id other) const;
  Id operator-(uword other) const;
  Id operator-(Id other) const;

  bool operator==(const Id &rhs) const;
  bool operator!=(const Id &rhs) const;

  bool operator<(const Id &rhs) const;
  bool operator>(const Id &rhs) const;
  bool operator<=(const Id &rhs) const;
  bool operator>=(const Id &rhs) const;

  friend std::ostream &operator<<(std::ostream &os, const Id &id);

  [[nodiscard]] std::array<uword, 3> to_xyz() const;

  [[nodiscard]] bool on_boundary() const;

  [[nodiscard]] bool on_corner() const;

  /// If this point is not on the boundary, will return `Id::Invalid`
  [[nodiscard]] Id boundary_inner() const;

  [[nodiscard]] std::array<Id, 6> adjacent() const;
};

#endif //MODELING_ID_H
