#define ARMA_USE_SUPERLU 1

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "Id.h"

//#include <TPolyMarker3D.h>
//#include <TCanvas.h>
//#include <TEnv.h>
//#include <TROOT.h>
//#include <TApplication.h>
//#include <TRootCanvas.h>

using namespace arma;

// todo command line args to config ROOT?
bool ROOT_MODE = true;

void solve_voltages() {//  std::ifstream geom_file("../geometry_small_text.ovf");
  std::ifstream geom_file("../geom/geometry_magnetic_network_L40_r0.5.ovf");
  if (!geom_file.is_open()) {
    std::cout << "Failed to read geometry file" << std::endl;
    return;
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
  const uword size = is_conductor.size();
  std::cout << "is_conductor.size() = " << size << std::endl;

  // todo this will be from geometry file somehow
// map from point to voltage
  std::map<Id, double> applied_voltages{
      {Id(20, 20, 0),  -5},
      {Id(20, 20, 40), 5}
  };
//  std::map<Id, double> applied_voltages{
//      {Id(1, 5, 0),  -5},
//      {Id(1, 5, SizeZ - 1), 5}
//  };

  if (ROOT_MODE) {
//  argc = 0;
//  argv = nullptr;
//  auto *app = new TApplication("app", &argc, argv);
//  auto *c = new TCanvas("c", "c");

//    auto *conductor_points = new TPolyMarker3D();
//    for (int id = 0; id < size; ++id) {
//      if (is_conductor[id]) {
//        auto xyz = Id(id).to_xyz();
//        conductor_points->SetNextPoint(xyz[0], xyz[1], xyz[2]);
//      }
//    }
//    std::cout << "conductor_points.Size = " << conductor_points->Size() << std::endl;
//    conductor_points->SetMarkerStyle(kFullDotLarge);
//    conductor_points->Draw();

//  c->Show();

//  c->Modified();
//  c->Update();

//  app->Run();
  }

  sp_mat A(size, size);
  vec b(size);
  for (uword i = 0; i < size; ++i) {
    const Id id(i);
    const auto idxyz = id.to_xyz();
    if (applied_voltages.find(id) != applied_voltages.end()) {
      // applied voltages are conductors (since voltage is constant)
      A(i, i) = 1;
      b(i) = applied_voltages[id];
    } else if (is_conductor[i] || id.on_corner()) {
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

//  A.brief_print("A:");

//  mat(A).print("A:");

  // todo: ±5V an/cathode, plot V -> E -> J

  superlu_opts opts;
  opts.symmetric = true;
  // todo calibrate this
  opts.pivot_thresh = 0.5;

  mat X = spsolve(A, b, "superlu", opts);
//  X.print("X: ");
  X.save("../data/top_bottom_voltages_L40_r0.5.bin");
  X.save(csv_name("../data/top_bottom_voltages_L40_r0.5.csv"));
}

void solve_e_field() {
  // for some reason loading this as a cube doesn't work
  mat voltage_data;
  voltage_data.load("../data/top_bottom_voltages_L40_r0.5.bin");

  cube voltages(SizeX, SizeY, SizeZ);
  {
    int z = 0;
    voltages.each_slice([&](mat &slice) {
      mat plane(mat(voltage_data(span(z, z + SizeX * SizeY), 0)));
      plane.reshape(SizeX, SizeY);
      slice = plane;
    });
  }

//  cube Ex(SizeX - 1, SizeY - 1, SizeZ - 1);
//  cube Ey(SizeX - 1, SizeY - 1, SizeZ - 1);
//  cube Ez(SizeX - 1, SizeY - 1, SizeZ - 1);


  for (uword z = 0; z < SizeZ; ++z) {
    mat Vx0 = voltages.slice(z)(span(0, SizeX - 1), span::all);
    mat Vx1 = voltages.slice(z)(span(1, SizeX), span::all);

    mat Ex = Vx1 - Vx1;
    Ex.print("Ex");
    return;
  }

//  mat col
//  voltages.cols(1, 1);
}

int main() {
  solve_voltages();
//  solve_e_field();

//  mat voltages;
//  voltages.load("../data/top_bottom_voltages_L40_r0.5.bin");
//  std::cout << "binCube.size() = " << voltages.size() << std::endl;
//  std::cout << "binCube.has_nan() = " << voltages.has_nan() << std::endl;
//
//  cube X(41, 41, 41);
//  {
//    uword i = 0;
////    X.each_slice([](mat &slice) { slice = })
//    mat eT('j',)
//  }

//  binCube.load("../data/top_bottom_voltages_L40_r0.5.bin");
//  std::cout << "binCube.size() = " << binCube.size() << std::endl;
//  std::cout << "binCube.has_nan() = " << binCube.has_nan() << std::endl;
}
