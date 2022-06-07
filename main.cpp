#include <iostream>
#include <array>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>

#include <TPolyMarker3D.h>
#include <thread>
#include <chrono>
#include <TCanvas.h>
#include <TEnv.h>

using namespace arma;

constexpr uword SizeX = 41;
constexpr uword SizeY = 41;
constexpr uword SizeZ = 41;
constexpr uword Sizes[3]{SizeX, SizeY, SizeZ};

struct Id {
  uword id;

  explicit Id(uword id) : id(id) {}

  Id(uword x, uword y, uword z) : id(x + y * SizeX + z * SizeX * SizeY) {}

  static Id Invalid;

  [[nodiscard]] bool is_valid() const {
    return id != Invalid.id;
  }

  Id operator+(uword other) const {
    return Id(id + other);
  }

  Id operator+(Id other) const {
    return *this + other.id;
  }

  Id operator-(uword other) const {
    return Id(id - other);
  }

  Id operator-(Id other) const {
    return *this - other.id;
  }

  friend std::ostream &operator<<(std::ostream &os, const Id &id) {
    const auto xyz = id.to_xyz();
    os << "Id(" << id.id << "), " << "(" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")";
    return os;
  }

  [[nodiscard]]
  std::array<uword, 3> to_xyz() const {
    std::array<uword, 3> ret = {this->id % SizeX, (this->id / SizeX) % SizeY, (this->id / SizeX / SizeY) % SizeZ};
    return ret;
  }

  [[nodiscard]]
  bool on_boundary() const {
    auto xyz = to_xyz();
    for (int i = 0; i < xyz.size(); ++i) {
      uword coord = xyz[i];
      if (coord == 0 || coord == Sizes[i] - 1)
        return true;
    }
    return false;
  }

  /// If this point is not on the boundary, will return `Id::Invalid`
  [[nodiscard]]
  Id boundary_inner() const {
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

  [[nodiscard]]
  std::array<Id, 6> adjacent() const {
    return {*this + 1, *this - 1, *this + SizeX, *this - SizeX, *this + SizeX * SizeY, *this - SizeX * SizeY};
  }
};

Id Id::Invalid = Id(std::numeric_limits<uword>::max());

//struct Point {
//  Id id;
//};

// todo command line args to config ROOT?
bool ROOT_MODE = false;

int main() {
  std::ifstream geom_file("../geometry_magnetic_network_L40_r0.5.ovf");
  if (!geom_file.is_open()) {
    std::cout << "Failed to read geometry file" << std::endl;
    return 1;
  }

  std::string line("#");
  for (; line.find('#') == 0; std::getline(geom_file, line)) {}

  std::vector<bool> is_conductor;
  unsigned int nconducting = 0;

  for (; line.find('#') != 0; std::getline(geom_file, line)) {
    bool conducts = line.find('0') != 0;
    is_conductor.push_back(conducts);
    nconducting += conducts;
  }
  std::cout << "is_conductor.size() = " << is_conductor.size() << std::endl;

  if (ROOT_MODE) {
    auto *conductor_points = new TPolyMarker3D();
    for (int id = 0; id < is_conductor.size(); ++id) {
      if (is_conductor[id]) {
        auto xyz = Id(id).to_xyz();
        conductor_points->SetNextPoint(xyz[0], xyz[1], xyz[2]);
      }
    }
    std::cout << "conductor_points.Size = " << conductor_points->Size() << std::endl;
    conductor_points->SetMarkerStyle(kFullDotLarge);
    conductor_points->Draw();
  }

  sp_mat A(is_conductor.size(), is_conductor.size());
  for (uword i = 0; i < is_conductor.size(); ++i) {
    const Id id(i);
    const auto idxyz = id.to_xyz();
    if (is_conductor[i]) {
      A(i, i) = 1;
    } else if (id.on_boundary()) {
      A(i, i) = 1;
      A(i, id.boundary_inner().id) = -1;
    } else {
      A(i, i) = -6;
      for (const auto &adjacent: id.adjacent()) {
        const auto adjxyz = adjacent.to_xyz();
        A(i, adjacent.id) = 1;
      }
    }
  }

  A.brief_print("A:");

//  mat X = spsolve(A, B);
}
