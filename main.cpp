#define ARMA_USE_SUPERLU 1

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "Id.h"
#include "utils.h"

//#include <TPolyMarker3D.h>
//#include <TCanvas.h>
//#include <TEnv.h>
//#include <TROOT.h>
//#include <TApplication.h>
//#include <TRootCanvas.h>

using namespace arma;

// todo command line args to config ROOT?
bool ROOT_MODE = true;
constexpr double h = 1.0;

struct Geometry {
  /// same size as `metal`
  std::vector<vec3> spins;
  /// each Id of a location in metal
  std::vector<Id> metal;
  /// if none, not metal. if not none, the index in `metal_ids` of this metal
  std::vector<std::optional<uword>> metal_idx;

  static Geometry from_file(const std::string &path) {
    std::ifstream geom_file(path);
    if (!geom_file.is_open()) {
      throw std::ios_base::failure(std::string("Failed to read geometry file '").append(path).append("'"));
    }

    std::string line("#");
    for (; line.find('#') == 0; std::getline(geom_file, line)) {}

    std::vector<vec3> spins;
    std::vector<Id> metal;
    std::vector<std::optional<uword>> metal_idx;

    for (; line.find('#') != 0; std::getline(geom_file, line)) {
      vec spin_growable(line);
      if (spin_growable.size() != 3) {
        spin_growable.print("spin");
        // todo throw better error
        throw "ERROR";
      }
      vec3 spin(spin_growable);
      if (!spin.is_zero()) {
        spins.push_back(spin);
        metal.emplace_back(metal_idx.size());
        metal_idx.emplace_back(metal.size() - 1);
      } else {
        metal_idx.emplace_back(std::nullopt);
      }
    }
    return {spins, metal, metal_idx};
  }
};

struct Solver {
  std::string name;
  Geometry geom;

  mat voltage;
  mat e_field;

  explicit Solver(std::string &&geometry_name)
      : geom(Geometry::from_file(std::string("../geom/").append(name).append(".ovf"))),
        name(std::move(geometry_name)) {}

  void solve_voltages() {
    const auto [_, metal, metal_idx] = geom;

    const uword size = metal.size();
    std::cout << "metal.size() = " << size << std::endl;

    // todo this will be from geometry file somehow
// map from point to voltage
    std::map<Id, double> applied_voltages;
    for (int y = 0; y < SizeY; ++y) {
      for (int x = 0; x < SizeX; ++x) {
      if (x == 1 && y == 1) continue;
        applied_voltages.emplace(Id(x, y, 0), 1);
        applied_voltages.emplace(Id(x, y, SizeZ - 1), 0);
      }
    }

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
            if (other && metal_idx[other->id]) {
              A(i, *metal_idx[other->id]) = 1;
              adj++;
            }
          }
        }
        // -6 from Poisson eq, then +1 for each Neumann boundary, for a total of -n_adj_metal
        A(i, i) = -adj;
      }
    }

    superlu_opts opts;
    opts.symmetric = true;
    // todo calibrate this
    opts.pivot_thresh = 0.5;

    voltage = spsolve(A, b, "superlu", opts);
    auto file = std::string("../data/voltages_").append(name);
    voltage.save(std::string(file).append(".bin"));
    voltage.save(csv_name(file.append(".csv")));
  }

  void solve_e_field() {
    const auto &[_, metal, metal_idx_clion_dumb] = geom;
    // for some reason CLion doesn't understand that lambda's can refer to destructing bindings
    const std::vector<std::optional<uword>> metal_idx = metal_idx_clion_dumb;

    const auto V = [&](Id id, Id backup) {
      return id.id > metal_idx.size() || !metal_idx[id.id]
             ? voltage(*metal_idx[backup.id])
             : voltage(*metal_idx[id.id]);
    };

    e_field = mat(metal.size(), 3);
    for (int i = 0; i < metal.size(); ++i) {
      const Id &id = metal[i];
      const auto [x_minus, x_plus, y_minus, y_plus, z_minus, z_plus] = id.adjacent();

      e_field(i, span::all) = {
          (V(x_plus, id) - V(x_minus, id)) / (2 * h),
          (V(y_plus, id) - V(y_minus, id)) / (2 * h),
          (V(z_plus, id) - V(z_minus, id)) / (2 * h)
      };
    }

    e_field.print("E");
    e_field.save(csv_name(std::string("../data/efield_").append(name).append(".csv")));
  }

  void solve_current() {
    const auto &[spin, metal, metal_idx_clion_dumb] = geom;

  }
};

int main() {
  Solver solver("wires_weird");
  solver.solve_voltages();
  solver.solve_e_field();
  solver.solve_current();
}
