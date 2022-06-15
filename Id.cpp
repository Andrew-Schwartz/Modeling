#include <armadillo>
#include <chrono>
#include <string>
#include <array>
#include <iostream>
#include "Id.h"

Id Id::Invalid = Id(std::numeric_limits<arma::uword>::max());

Id::Id(uword id) : id(id) {}

Id::Id(uword x, uword y, uword z) : id(x + y * SizeX + z * SizeX * SizeY) {}

Id::Id(const uvec3 &xyz) : Id(xyz[0], xyz[1], xyz[2]) {}

bool Id::is_valid() const {
  return id != Invalid.id;
}

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
    if (std::signbit(disp) && std::abs(disp) > coord || !std::signbit(disp) && coord + disp >= Sizes[i]) {
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
  uvec3 ret = {this->id % SizeX, (this->id / SizeX) % SizeY, (this->id / SizeX / SizeY) % SizeZ};
  return ret;
}

bool Id::on_boundary() const {
  auto xyz = to_xyz();
  for (int i = 0; i < xyz.size(); ++i) {
    uword coord = xyz[i];
    if (coord == 0 || coord == Sizes[i] - 1)
      return true;
  }
  return false;
}

bool Id::on_corner() const {
  const auto xyz = to_xyz();
  int boundaries = 0;
  for (int i = 0; i < xyz.size(); ++i) {
    uword coord = xyz[i];
    if (coord == 0 || coord == Sizes[i] - 1)
      boundaries++;
  }
  return boundaries >= 2;
}

Id Id::boundary_inner() const {
  auto xyz = to_xyz();
  uword x = xyz[0], y = xyz[1], z = xyz[2];
  Id ret = x == 0 ? *this + 1
                  : x == SizeX - 1 ? *this - 1
                                   : y == 0 ? *this + SizeX
                                            : y == SizeY - 1 ? *this - SizeX
                                                             : z == 0 ? *this + SizeX * SizeY
                                                                      : z == SizeZ - 1 ? *this - SizeX * SizeY
                                                                                       : Id::Invalid;
  if (!ret.is_valid()) {
    std::cout << "No boundary for " << this << std::endl;
    // todo
//      throw
  }
  return ret;
}

std::array<Id, 6> Id::adjacent() const {
  return {*this + 1, *this - 1, *this + SizeX, *this - SizeX, *this + SizeX * SizeY, *this - SizeX * SizeY};
}
