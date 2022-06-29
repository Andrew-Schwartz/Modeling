#include <armadillo>
#include <string>
#include <iostream>
#include "Id.h"

uword Id::SizeX = 0;
uword Id::SizeY = 0;
uword Id::SizeZ = 0;

Id::Id(uword id) : id(id) {}

Id::Id(uword x, uword y, uword z) : id(x + y * SizeX + z * SizeX * SizeY) {}

Id::Id(const uvec3 &xyz) : Id(xyz[0], xyz[1], xyz[2]) {}

Id Id::operator+(uword other) const {
  return Id(id + other);
}

Id Id::operator+(Id other) const {
  return *this + other.id;
}

std::optional<Id> Id::operator+(const ivec3 &displacement) const {
  uvec3 xyz = to_xyz();
  for (int i = 0; i < 3; ++i) {
    sword disp = displacement(i);
    uword &coord = xyz(i);
    if (std::signbit(disp) && std::abs(disp) > coord || !std::signbit(disp) && coord + disp >= Sizes()[i]) {
      return std::nullopt;
    }
    coord += disp;
  }
  return std::make_optional(Id(xyz));
}

Id Id::operator-(uword other) const {
  return Id(id - other);
}

Id Id::operator-(Id other) const {
  return *this - other.id;
}

bool Id::operator==(const Id &rhs) const {
  return id == rhs.id;
}

bool Id::operator!=(const Id &rhs) const {
  return !(rhs == *this);
}

bool Id::operator<(const Id &rhs) const {
  return id < rhs.id;
}

bool Id::operator>(const Id &rhs) const {
  return rhs < *this;
}

bool Id::operator<=(const Id &rhs) const {
  return !(rhs < *this);
}

bool Id::operator>=(const Id &rhs) const {
  return !(*this < rhs);
}

std::ostream &operator<<(std::ostream &os, const Id &id) {
  const auto xyz = id.to_xyz();
  os << "Id(" << id.id << "), " << "(" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")";
  return os;
}

uvec3 Id::to_xyz() const {
  uvec3 ret = {id % SizeX, (id / SizeX) % SizeY, (id / SizeX / SizeY) % SizeZ};
  return ret;
}

double Id::dist3d(const Id &other) const {
  uvec3 xyz = to_xyz(),
      oxyz = other.to_xyz();
  double sum_of_squares = 0;
  for (int i = 0; i < 3; ++i)
    sum_of_squares += std::pow(sword(xyz(i)) - sword(oxyz(i)), 2);
  return std::sqrt(sum_of_squares);
}