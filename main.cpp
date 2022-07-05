#define ARMA_USE_SUPERLU 1

#include <iostream>
#include <iomanip>
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

struct SolverOpts {
  /// ignore spin texture in input file and treat this direction as the unit vector spin for all magnetic points
  std::optional<vec3> uniform_spin;

  /// if one of the points in the voltage file is not in the network, fall back to the closest point(s) to it
  bool use_closest_point;

  /// uniform spin pointing from -s_hat to s_hat
  bool vary_spins;
};

class Solver {
  static std::string make_name(
      const std::string &geom_name,
      const std::string &voltage_name,
      const SolverOpts &opts
  ) {
    std::ostringstream name;
    name << geom_name << "_" << voltage_name;
//    auto name = std::string(geom_name).append("_").append(voltage_name);
    if (opts.uniform_spin) {
      name << std::fixed << std::setprecision(1);
      name << "_spin" << opts.uniform_spin->at(0)
           << "," << opts.uniform_spin->at(1)
           << "," << opts.uniform_spin->at(2);
    }
    return name.str();
  }

public:
  std::string name;
  Geometry geom;
  std::map<Id, double> applied_volts;
  bool save_files = true;
  vec3 s_hat;
  mat33 rotation;

  mat voltage;
  mat e_field;
  mat current;

  Solver(const std::string &geometry_name, const std::string &voltage_name, const SolverOpts &opts = SolverOpts())
      : name(make_name(geometry_name, voltage_name, opts)),
        geom(Geometry::from_file(std::string("../geom/").append(geometry_name).append(".ovf"))) {
    if (opts.uniform_spin) {
      std::fill(geom.spins.begin(), geom.spins.end(), *opts.uniform_spin);
    }

    auto path = std::string("../geom/volts_").append(voltage_name).append(".tsv");
    std::ifstream volt_file(path);
    if (!volt_file.is_open()) {
      throw std::ios_base::failure(std::string("Failed to read volt file '").append(path).append("'"));
    }

    // todo subclass exception
    std::vector<std::string> split;
    auto stoull = [&split](size_t idx) { return std::stoull(split[idx]); };
    auto stod = [&split](size_t idx) { return std::stod(split[idx]); };
    for (std::string line; std::getline(volt_file, line);) {
      split.clear();
      std::istringstream stream(line);
      std::string tmp;
      while (stream >> tmp) {
        split.push_back(tmp);
      }
      // id, volts
      if (split.size() == 2) {
        applied_volts.emplace(Id(stoull(0)), stod(1));
      } else if (split.size() == 4) {
        // todo finish this maybe if it will be useful
//        std::array<std::vector<uword>, 3> xyz;
//        for (size_t dim: {0, 1, 2}) {
//
//          auto str = split[dim];
//          if (str.find(">=") != std::string::npos) {
//            uword min = std::stoull(str.substr(std::strlen(">=")));
//            for (uword i = min; i < Id::Sizes()[dim]; ++i)
//              xyz[dim].push_back(i);
//          } else if (str.find(">") != std::string::npos) {
//            uword min = std::stoull(str.substr(std::strlen(">")));
//            for (uword i = min; i < Id::Sizes()[dim]; ++i)
//              xyz[dim].push_back(i);
//          } else if (str.find("<=") != std::string::npos) {
//            uword min = std::stoull(str.substr(std::strlen("<=")));
//            for (uword i = 0; i < Id::Sizes()[dim]; ++i)
//              xyz[dim].push_back(i);
//          }
//        }
        applied_volts.emplace(Id(stoull(0), stoull(1), stoull(2)), stod(3));
      } else {
        throw std::runtime_error(std::string("Unable to parse voltage line: ").append(line));
      }
    }
    // bad id -> good id
    std::map<Id, Id> replacements;
    bool any_not_in_metal = false;
    for (const auto &[id, v]: applied_volts) {
      // if not in the metal, print the closest ids
      if (!geom.metal_idx[id.id]) {
        any_not_in_metal = true;
        std::cout << id << " is not a point in the network" << std::endl;
        double min_dist = std::numeric_limits<double>::max();
        std::vector<Id> closest;
        for (const Id &other: geom.metal) {
          const double dist = id.dist3d(other);
          if (dist < min_dist) {
            min_dist = dist;
            closest = {other};
          } else if (dist == min_dist) {
            closest.push_back(other);
          }
        }
        std::cout << (closest.size() > 1 ? "These are" : "This is")
                  << " the closest point"
                  << (closest.size() > 1 ? "s" : "")
                  << ", "
                  << min_dist
                  << " units away:"
                  << std::endl;
        for (const Id &close: closest) {
          std::cout << close << std::endl;
        }
        if (opts.use_closest_point) {
          Id use = closest[0];
          std::cout << (closest.size() > 1 ? "The first" : "This") << " point will be used instead\n" << std::endl;
          replacements.emplace(id, use);
        } else {
          std::cout << std::endl;
        }
      }
    }
    if (opts.use_closest_point) {
      for (const auto &[bad, good]: replacements) {
        auto node = applied_volts.extract(bad);
        node.key() = good;
        applied_volts.insert(std::move(node));
      }
    } else if (any_not_in_metal && !opts.use_closest_point) {
      throw std::runtime_error(
          "Invalid points used in voltage configuration. Choose different points or enable `SolverOpts.use_closest_point`."
      );
    }
    if (opts.vary_spins) {
      save_files = false;
      if (applied_volts.size() != 2
//      || std::none_of(applied_volts.begin(), applied_volts.end(), [](const auto &[id, v]) { return v == 0; })
          ) {
        throw std::runtime_error("Can only vary spins with 1 input and 1 output port");
      }
      Id input(0), output(0);
      for (const auto &[id, v]: applied_volts) {
        if (v == 0) output = id;
        else input = id;
      }
      vec3 s = conv_fixed<3, uword, double>(output.to_xyz() - input.to_xyz());
      s_hat = normalise(s);
      std::cout << "s_hat = " << s_hat << std::endl;
    }
  }

  void save(const mat &data, const csv_name &file_name) const {
    if (!save_files) return;
    const mat &save = join_rows(geom.xyz(), data);
    save.save(file_name);
  }

  void solve_voltages() {
    const auto [_, metal, metal_idx] = geom;

    const uword size = metal.size();
    std::cout << "metal.size() = " << size << std::endl;

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

    field<std::string> header{"x", "y", "z", "jx", "jy", "jz"};
    save(current, csv_name(std::string("../data/current_").append(name).append(".csv"), header));
  }

  void rotate_spin() {
    // perpendicular to s_hat so we can rotate around it
    vec3 rotation_axis {s_hat(1), -s_hat(0), 0};
    rotation_axis = normalise(rotation_axis);

    Id output_id = std::find_if(applied_volts.cbegin(), applied_volts.cend(), [](const auto &pair) { return pair.second == 0; })->first;
    uword output_idx = *geom.metal_idx[output_id.id];

    const int n_rotations = 360;
    mat current_magnitude(n_rotations, 2); 
    for (int degrees = 0; degrees < n_rotations; ++degrees) {
      const double theta = degrees * datum::pi / 180.0;
      vec3 spin = rotate_vector(rotation_axis, theta, s_hat);
      std::fill(geom.spins.begin(), geom.spins.end(), spin);
      solve_current();
      rowvec out = current(output_idx, span::all);
      current_magnitude(degrees, span::all) = {theta, norm(out)};
    }
    field<std::string> header{"theta", "|j|"};
    current_magnitude.save(csv_name(join_string("../data/current_varied_", name, ".csv"), header));
  }
};

// 1 in 1 out
// vector n from in to out, rotate spin from negative n to positive n

int main() {
  Solver solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners",
                SolverOpts{.use_closest_point = true, .vary_spins = true});
  solver.solve_voltages();
  solver.solve_e_field();

  solver.rotate_spin();

  /*
   * 0.0: 0.0055   0.0059   0.0079
   * 0.1: 0.0102   0.0097   0.0153
   *   z: 0.0061   0.0063   0.0132
   *   y: 0.0061   0.0126   0.0066
   *   x: 0.0123   0.0063   0.0066
   */
  /*
  std::vector<Solver> solvers{
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners", SolverOpts{std::make_optional(vec3{0, 0, 1}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners", SolverOpts{std::make_optional(vec3{0, 0, -1}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners", SolverOpts{std::make_optional(vec3{0, 1, 0}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners", SolverOpts{std::make_optional(vec3{1, 0, 0}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_corners", SolverOpts{.use_closest_point = true}),
      Solver("spin_wire4current_1_B=0.1", "sw4c10.0_corners", SolverOpts{.use_closest_point = true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_three", SolverOpts{std::make_optional(vec3{0, 0, 1}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_three", SolverOpts{std::make_optional(vec3{0, 0, -1}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_three", SolverOpts{std::make_optional(vec3{0, 1, 0}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_three", SolverOpts{std::make_optional(vec3{1, 0, 0}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_three", SolverOpts{.use_closest_point = true}),
      Solver("spin_wire4current_1_B=0.1", "sw4c10.0_three", SolverOpts{.use_closest_point = true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_more", SolverOpts{std::make_optional(vec3{0, 0, 1}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_more", SolverOpts{std::make_optional(vec3{0, 0, -1}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_more", SolverOpts{std::make_optional(vec3{0, 1, 0}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_more", SolverOpts{std::make_optional(vec3{1, 0, 0}), true}),
      Solver("spin_wire4current_1_B=0.0", "sw4c10.0_more", SolverOpts{.use_closest_point = true}),
      Solver("spin_wire4current_1_B=0.1", "sw4c10.0_more", SolverOpts{.use_closest_point = true}),
  };

  std::ofstream out("../data/spin_wires4current_1_output_current.txt");

  for (Solver solver: std::move(solvers)) {
    solver.solve_voltages();
    solver.solve_e_field();
    solver.solve_current();
    out << "-----------------------\n"
        << solver.name
        << std::endl;
    for (const auto &[id, volts]: solver.applied_volts) {
      if (volts == 0) {
        auto current_row = solver.current(*solver.geom.metal_idx[id.id], span::all);
        out << id << ", current = " << current_row;
      }
    }
  }
  out << "-----------------------" << std::endl;*/

/*  SolverOpts opts{
//    .uniform_spin = std::make_optional(vec3 {0, 0, 1})
//    .uniform_spin = std::make_optional(vec3 {0, 1, 0})
//    .uniform_spin = std::make_optional(vec3 {1, 0, 0})
      .use_closest_point = true
  };

//  Solver solver("wires", "four_wires_all");
//  Solver solver("cube", "cube555corners");
//  Solver solver("1d", "1d", opts);
  Solver solver("spin_wire4current_3_B=0.0", "sw4c10.0_corners", opts);
  solver.solve_voltages();
  solver.solve_e_field();
  solver.solve_current();

  for (const auto &[id, volts]: solver.applied_volts) {
    if (volts == 0) {
      auto current_row = solver.current(*solver.geom.metal_idx[id.id], span::all);
      std::cout << id << ", current = " << current_row;
    }
  }*/
}
