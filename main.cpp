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

void solve_voltages() {
  std::ifstream geom_file("../geom/wires.ovf");
//  std::ifstream geom_file("../geom/geometry_magnetic_network_L40_r0.5.ovf");
  if (!geom_file.is_open()) {
    std::cout << "Failed to read geometry file" << std::endl;
    return;
  }

  std::string line("#");
  for (; line.find('#') == 0; std::getline(geom_file, line)) {}

  std::vector<bool> is_conductor;
  std::vector<Id> metal;

  for (; line.find('#') != 0; std::getline(geom_file, line)) {
    bool conducts = line.find('0') != 0;
    if (conducts) {
      metal.emplace_back(is_conductor.size());
    }
    is_conductor.push_back(conducts);
  }
  const uword size = metal.size();
  std::cout << "metal.size() = " << size << std::endl;

  // todo this will be from geometry file somehow
// map from point to voltage
  std::map<Id, double> applied_voltages{
//      {Id(1, 1, 0),  0},
//      {Id(1, 1, 2), 1}
  };
  for (int y = 0; y < SizeY; ++y) {
    for (int x = 0; x < SizeX; ++x) {
      applied_voltages.emplace(Id(x, y, 0), 1);
      applied_voltages.emplace(Id(x, y, SizeZ - 1), 0);
    }
  }
//  std::map<Id, double> applied_voltages{
//      {Id(1, 5, 0),  -5},
//      {Id(1, 5, SizeZ - 1), 5}
//  };

//  if (ROOT_MODE) {
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
//  }

  sp_mat A(size, size);
  vec b(size);
  for (uword i = 0; i < size; ++i) {
    const Id id = metal[i];
    if (applied_voltages.find(id) != applied_voltages.end()) {
      A(i, i) = 1;
      b(i) = applied_voltages[id];
    } else {
      int adj = 0;
      for (const sword offset: {-1, 1}) {
        for (const uword dim: {0, 1, 2}) {
          ivec3 displacement;
          displacement(dim) = offset;
          auto other = id + displacement;
          if (other) {
            if (is_conductor[other->id]) {
              uword metal_idx = std::find(metal.begin(), metal.end(), other) - metal.begin();
              A(i, metal_idx) = 1;
              adj++;
            }
          }
        }
      }
      // because some points aren't in the domain, its not always -6 here (todo prove mathematically, but this *is*
      //  correct because otherwise voltage doesn't propagate properly, due to it trying to go in all 6 directions)
      A(i, i) = -adj;
    }
    /* else if (is_conductor[i] || id.on_corner()) {
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
    }*/
  }

//  A.brief_print("A:");

  mat(A).print("A:");
  b.print("b: ");

  // todo: ±5V an/cathode, plot V -> E -> J

  superlu_opts opts;
  opts.symmetric = true;
  // todo calibrate this
  opts.pivot_thresh = 0.5;

  mat X = spsolve(A, b, "superlu", opts);
//  X.print("X: ");
  std::string file("../data/voltages_wires");
  X.save(std::string(file).append(".bin"));
  X.save(csv_name(file.append(".csv")));
}

void solve_e_field() {
  // for some reason loading this as a cube doesn't work
  mat voltage_data;
  voltage_data.load("../data/top_bottom_voltages_L40_r0.5.bin");

  cube voltages(SizeX, SizeY, SizeZ);
  {
    int z = 0;
    voltages.each_slice([&](mat &slice) {
      mat plane(mat(voltage_data(span(z, z
      +SizeX * SizeY), 0)));
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

    mat Ex = Vx1 - Vx0;
    Ex.print("Ex");
    return;
  }

//  mat col
//  voltages.cols(1, 1);
}

int main() {
  solve_voltages();
//  solve_e_field();
}
