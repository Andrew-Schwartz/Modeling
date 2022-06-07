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

using namespace arma;
using usize = unsigned int;

constexpr usize SizeX = 41;
constexpr usize SizeY = 41;
constexpr usize SizeZ = 41;

struct Id {
  usize id;

  explicit Id(usize id) : id(id) {}

  Id(usize x, usize y, usize z) : id(x + y * (SizeX + 1) + z * (SizeX + 1) * (SizeY + 1)) {}

  [[nodiscard]]
  std::array<usize, 3> to_xyz() const {
    constexpr usize sx = SizeX + 1;
    constexpr usize sy = SizeY + 1;
    constexpr usize sz = SizeZ + 1;
    std::array<usize, 3> ret = {this->id % sx, (this->id / sx) % sy, (this->id / sx / sy) % sz};
    return ret;
  }
};

//struct Point {
//  Id id;
//};

int main() {
  auto *c1 = new TCanvas("c1", "c1");

  std::ifstream geom_file("geometry_magnetic_network_L40_r0.5.ovf");
  if (!geom_file.is_open()) {
    std::cout << "Failed to read geometry file" << std::endl;
    return 1;
  }

  std::string line("#");
  for (; line.find('#') == 0; std::getline(geom_file, line)) {}

  std::vector<bool> conductor;

  for (; line.find('#') != 0; std::getline(geom_file, line)) {
    conductor.push_back(line.find('0') != 0);
  }
  std::cout << "conductor.size() = " << conductor.size() << std::endl;

//  arma::sp_mat

  auto *conductor_points = new TPolyMarker3D();
  for (int id = 0; id < conductor.size(); ++id) {
    if (conductor[id]) {
      auto xyz = Id(id).to_xyz();
      conductor_points->SetNextPoint(xyz[0], xyz[1], xyz[2]);
    }
  }
  std::cout << "conductor_points.Size = " << conductor_points->Size() << std::endl;
  conductor_points->SetMarkerStyle(kFullDotLarge);
  conductor_points->Draw();
//  c.Draw();
//  c.Show();
}
