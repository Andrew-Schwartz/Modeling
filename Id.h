#ifndef MODELING_ID_H
#define MODELING_ID_H

#include <armadillo>

using namespace arma;

constexpr uword SizeX = 5; // 41; // 20;
constexpr uword SizeY = 5; // 41; // 20;
constexpr uword SizeZ = 10; // 41; // 40;
constexpr uword Sizes[3]{SizeX, SizeY, SizeZ};

struct Id {
  uword id;

  explicit Id(uword id);
  Id(uword x, uword y, uword z);
  explicit Id(const uvec3 &xyz);

  static Id Invalid;

  [[nodiscard]] bool is_valid() const;

  Id operator+(uword other) const;
  Id operator+(Id other) const;
  /// returns std::nullopt if this operation would result in an Id outside of the sample
  std::optional<Id> operator+(const ivec3& displacement) const;
  Id operator-(uword other) const;
  Id operator-(Id other) const;

  bool operator==(const Id &rhs) const;
  bool operator!=(const Id &rhs) const;

  bool operator<(const Id &rhs) const;
  bool operator>(const Id &rhs) const;
  bool operator<=(const Id &rhs) const;
  bool operator>=(const Id &rhs) const;

  friend std::ostream &operator<<(std::ostream &os, const Id &id);

  [[nodiscard]] uvec3 to_xyz() const;

  [[nodiscard]] bool on_boundary() const;

  [[nodiscard]] bool on_corner() const;

  /// If this point is not on the boundary, will return `Id::Invalid`
  [[nodiscard]] Id boundary_inner() const;

  /// x, y, z order with -1 before +1
  [[nodiscard]] std::array<Id, 6> adjacent() const;
};

#endif //MODELING_ID_H
