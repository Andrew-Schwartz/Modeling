#define ARMA_USE_SUPERLU 1

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include "Id.h"
#include "utils.h"

using namespace arma;

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
    for (; line.find('#') == 0; std::getline(geom_file, line)) {
      if (line.find("nodes:") != std::string::npos) {
        uword size = std::stoll(line.substr(line.rfind(' ')));
        if (line.find('x') != std::string::npos)
          Id::SizeX = size;
        else if (line.find('y') != std::string::npos)
          Id::SizeY = size;
        else if (line.find('z') != std::string::npos)
          Id::SizeZ = size;
      }
    }

    if (Id::SizeX == 0 || Id::SizeY == 0 || Id::SizeZ == 0) {
      throw std::runtime_error("Did not set size of sample");
    }

    std::vector<vec3> spins;
    std::vector<Id> metal;
    std::vector<std::optional<uword>> metal_idx;

    for (; line.find('#') != 0; std::getline(geom_file, line)) {
      vec spin_growable(line);
      if (spin_growable.size() != 3) {
        spin_growable.print("spin");
        // todo throw better error
        throw std::runtime_error("ERROR");
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

    auto sizes = Id::Sizes();
    if (metal_idx.size() != std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<>())) {
      throw std::runtime_error("Id dimensions are wrong!");
    }

    return {spins, metal, metal_idx};
  }

  /// x, y, z coordinates in order by metal id
  [[nodiscard]] mat xyz() const {
    mat coords(metal.size(), 3);
    for (int i = 0; i < metal.size(); ++i) {
      umat xyz = umat(Id(metal[i]).to_xyz());
      xyz.reshape(1, 3);
      coords(i, span::all) = conv_to<mat>::from(xyz);
    }
    return coords;
  }

  [[nodiscard]] std::array<std::optional<Id>, 6> adjacent(const Id &id) const {
    std::array<std::optional<Id>, 6> ret;
    int i = 0;
    for (const uword dim: {0, 1, 2}) {
      for (const sword offset: {-1, 1}) {
        ivec3 displacement;
        displacement(dim) = offset;
        auto adj = id + displacement;
        if (adj && metal_idx[adj->id]) {
          ret[i] = adj;
        }
        ++i;
      }
    }
    return ret;
  }
};

struct Solver {
  std::string name;
  Geometry geom;
  std::map<Id, double> applied_volts;

  mat voltage;
  mat e_field;
  mat current;

  explicit Solver(const std::string &geometry_name, const std::string &voltage_name)
      : name(std::string(geometry_name).append("_").append(voltage_name)),
        geom(Geometry::from_file(std::string("../geom/").append(geometry_name).append(".ovf"))) {
    auto path = std::string("../geom/volts_").append(voltage_name).append(".tsv");
    std::ifstream volt_file(path);
    if (!volt_file.is_open()) {
      throw std::ios_base::failure(std::string("Failed to read volt file '").append(path).append("'"));
    }

    uword x, y, z;
    double volts;
    for (std::string line; std::getline(volt_file, line);) {
      std::istringstream(line) >> x >> y >> z >> volts;
      applied_volts.emplace(Id(x, y, z), volts);
    }
  }

  void save(const mat &data, const csv_name &name) const {
    const mat &save = join_rows(geom.xyz(), data);
    save.save(name);
  }

  void solve_voltages() {
    const auto [_, metal, metal_idx] = geom;

    const uword size = metal.size();
    std::cout << "metal.size() = " << size << std::endl;

    // todo this will be from geometry file somehow
// map from point to voltage
//    std::map<Id, double> applied_voltages{
//        {Id(0, 0, 0),                         1},
//        {Id(SizeX - 1, SizeY - 1, SizeZ - 1), 0},
//    };
//    applied_volts = applied_voltages;

    sp_mat A(size, size);
    vec b(size);
    for (uword i = 0; i < size; ++i) {
      const Id id = metal[i];
      if (applied_volts.find(id) != applied_volts.end()) {
        A(i, i) = 1;
        b(i) = applied_volts[id];
      } else {
        int n_adj = 0;
        for (const auto adj: geom.adjacent(id)) {
          if (adj) {
            A(i, *metal_idx[adj->id]) = 1;
            ++n_adj;
          }
        }
        // -6 from Poisson eq, then +1 for each Neumann boundary, for a total of -n_adj_metal
        A(i, i) = -n_adj;
      }
    }

//    mat(A).print("A");
//    b.print("b");

    superlu_opts opts;
    opts.symmetric = true;
    // todo calibrate this
    opts.pivot_thresh = 0.5;

    voltage = spsolve(A, b, "superlu", opts);
    field<std::string> header{"x", "y", "z", "volts"};
    save(voltage, csv_name(std::string("../data/voltages_").append(name).append(".csv"), header));
  }

  void solve_e_field() {
    const auto &[_, metal, metal_idx_clion_dumb] = geom;
    // for some reason CLion doesn't understand that lambda's can refer to destructing bindings
    const std::vector<std::optional<uword>> &metal_idx = metal_idx_clion_dumb;

    const auto V = [&](std::optional<Id> id, Id backup) {
      return id
             ? voltage(*metal_idx[id->id])
             : voltage(*metal_idx[backup.id]);
    };

    e_field = mat(metal.size(), 3);
    for (int i = 0; i < metal.size(); ++i) {
      const Id &id = metal[i];
      const auto [x_minus, x_plus, y_minus, y_plus, z_minus, z_plus] = geom.adjacent(id);

      e_field(i, span::all) = {
          (V(x_minus, id) - V(x_plus, id)) / (2 * h),
          (V(y_minus, id) - V(y_plus, id)) / (2 * h),
          (V(z_minus, id) - V(z_plus, id)) / (2 * h)
      };
    }

//    e_field.print("E");
//    e_field.save(csv_name(std::string("../data/efield_").append(name).append(".csv")));
    field<std::string> header{"x", "y", "z", "Ex", "Ey", "Ez"};
    save(e_field, csv_name(std::string("../data/efield_").append(name).append(".csv"), header));
  }

  void solve_current() {
    const auto &[spins, metal, metal_idx] = geom;

    double rho_parallel = 5;
    double rho_perpendicular = 10;
    double sigma_nought = (1.0 / rho_parallel + 2.0 / rho_perpendicular) / 3.0;
    double sigma_one = (1.0 / rho_parallel - 1.0 / rho_perpendicular);

    current = mat(metal.size(), 3);
    for (int i = 0; i < metal.size(); ++i) {
      const auto &spin = spins[i];
      double mx = spin(0), my = spin(1), mz = spin(2);
      // @formatter:off
      mat projection{
          {mx * mx - 1.0 / 3.0, mx * my            , mx * mz            },
          {my * mx            , my * my - 1.0 / 3.0, my * mz            },
          {mz * mx            , mz * my            , mz * mz - 1.0 / 3.0}
      };
      // @formatter:on
      const mat sigma = sigma_nought * eye(3, 3) + sigma_one * projection;

      mat j = sigma * reshape(e_field(i, span::all), 3, 1);
      j.reshape(1, 3);
      current(i, span::all) = j;
    }

    field<std::string> header{"x","y","z","jx","jy","jz"};
    save(current, csv_name(std::string("../data/current_").append(name).append(".csv"), header));
  }
};

int main() {
//  Solver solver("wires_weird2", "four_wires_three");
//  Solver solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners");
//  Solver solver("cube", "cube555corners");
  Solver solver("1d", "1d");
  solver.solve_voltages();
  solver.solve_e_field();
  solver.solve_current();
}
