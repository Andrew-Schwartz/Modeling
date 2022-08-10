#define ARMA_USE_SUPERLU 1

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <filesystem>
#include <armadillo>
#include "Id.h"
#include "utils.h"
#include "ovf_wrapper.hpp"

using namespace arma;

struct Geometry {
    double h;
    /// same size as `metal`
    std::vector<vec3> spins;
    /// each Id of a location in metal
    std::vector<Id> metal;
    /// if none, not metal. if not none, the index in `metal_ids` of this metal
    std::vector<std::optional<size_t>> metal_idx;

    static Geometry from_file(const std::string &path) {
        auto file = ovf::File::open_file(path.c_str());
        const mat data = file.read();
        auto &seg = file.segments[0];
        Id::SizeX = seg->n_cells[0];
        Id::SizeY = seg->n_cells[1];
        Id::SizeZ = seg->n_cells[2];
        double h = seg->step_size[0];

        std::vector<vec3> spins;
        std::vector<Id> metal;
        std::vector<std::optional<size_t>> metal_idx;

        data.each_row([&](const rowvec &row) {
            vec3 spin = row.t();
            if (!spin.is_zero()) {
                spins.push_back(spin);
                metal.emplace_back(metal_idx.size());
                metal_idx.emplace_back(metal.size() - 1);
            } else {
                metal_idx.emplace_back(std::nullopt);
            }
        });

/*        std::ifstream geom_file(path);
        if (!geom_file.is_open()) {
            throw std::ios_base::failure(join_string("Failed to read geometry file '", path, "'"));
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
        }*/

        /*for (; line.find('#') != 0; std::getline(geom_file, line)) {
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
        }*/

        auto sizes = Id::Sizes();
        if (metal_idx.size() != std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<>())) {
            throw std::runtime_error("Id dimensions are wrong!");
        }

        return {h, spins, metal, metal_idx};
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

    // for timing:
    double solve_linear_eq_timing = 0;
    double io_timing = 0;
    double solve_efield_timing = 0;

    mat voltage;
    mat e_field;
    mat current;

    Solver(const std::string &geometry_name, const std::string &voltage_name, const SolverOpts &opts = SolverOpts())
            : name(make_name(geometry_name, voltage_name, opts)),
              geom(Geometry::from_file(join_string("../geom/", geometry_name, ".ovf"))) {
        std::cout << "Solving " << name << std::endl;

        if (opts.uniform_spin) {
            std::fill(geom.spins.begin(), geom.spins.end(), *opts.uniform_spin);
        }

        auto path = join_string("../geom/volts_", voltage_name, ".tsv");
        std::ifstream volt_file(path);
        if (!volt_file.is_open()) {
            throw std::ios_base::failure(join_string("Failed to read volt file '", path, "'"));
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
                /*std::array<std::vector<uword>, 3> xyz;
                for (size_t dim: {0, 1, 2}) {

                  auto str = split[dim];
                  if (str.find(">=") != std::string::npos) {
                    uword min = std::stoull(str.substr(std::strlen(">=")));
                    for (uword i = min; i < Id::Sizes()[dim]; ++i)
                      xyz[dim].push_back(i);
                  } else if (str.find(">") != std::string::npos) {
                    uword min = std::stoull(str.substr(std::strlen(">")));
                    for (uword i = min; i < Id::Sizes()[dim]; ++i)
                      xyz[dim].push_back(i);
                  } else if (str.find("<=") != std::string::npos) {
                    uword min = std::stoull(str.substr(std::strlen("<=")));
                    for (uword i = 0; i < Id::Sizes()[dim]; ++i)
                      xyz[dim].push_back(i);
                  }
                }*/
                applied_volts.emplace(Id(stoull(0), stoull(1), stoull(2)), stod(3));
            } else {
                throw std::runtime_error(join_string("Unable to parse voltage line: ", line));
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
                    std::cout << (closest.size() > 1 ? "The first" : "This") << " point will be used instead\n"
                              << std::endl;
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
            Id input = std::find_if(applied_volts.begin(), applied_volts.end(),
                                    [](const auto &pair) { return pair.second != 0; })->first;
            Id output = std::find_if(applied_volts.begin(), applied_volts.end(),
                                     [](const auto &pair) { return pair.second == 0; })->first;

            vec3 s = conv_fixed<3, uword, double>(output.to_xyz() - input.to_xyz());
            s_hat = normalise(s);
            std::cout << "s_hat = " << s_hat << std::endl;
        }
    }

    void save(const mat &data, const csv_name &file_name) const {
        const mat save = join_rows(geom.xyz(), data);
        save.save(file_name);
    }

    void save_ovf(const mat &data, const char *filename) const {
        auto segment = ovf::Segment()
                .title("E")
                .mesh_type(ovf::MeshType::Rectangular)
                .mesh_unit("m")
                .bounds({0, 0, 0}, {1.55e-7, 1.55e-7, 1.55e-7})
                .value_dim(3)
                .value_labels("E_x E_y E_z")
                .value_units("V/m V/m V/m")
//                .point_count(geom.metal.size());
                .n_nodes({(int) Id::SizeX, (int) Id::SizeY, (int) Id::SizeZ})
                .base({2.5e-9, 2.5e-9, 2.5e-9})
                .step_size({5e-9, 5e-9, 5e-9});
        size_t size = geom.metal_idx.size();
        mat save(3, size);
        for (size_t i = 0; i < size; ++i) {
            auto idx = geom.metal_idx[i];
            save(span::all, i) = idx
                                 ? data(*idx, span::all).t()
                                 : vec{0, 0, 0};
        }
//        const mat save = join_rows(geom.xyz(), data).t();
        ovf::File::write_file(filename, std::move(segment), save.memptr());
    }

    void solve_voltages() {
        const auto [_1, _2, metal, metal_idx] = geom;

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
        solve_linear_eq_timing = timer::time([&]() { this->voltage = spsolve(A, b, "superlu", opts); });
        field<std::string> header{"x", "y", "z", "volts"};
        save(voltage, csv_name(join_string("../data/voltages_", name, ".csv"), header));
    }

    void solve_e_field() {
        const double h = geom.h;
        const auto &metal = geom.metal;
        const auto &metal_idx = geom.metal_idx;

        const auto V = [&](std::optional<Id> id, Id backup) {
            return 1e-6 * (id
                           ? voltage(*metal_idx[id->id])
                           : voltage(*metal_idx[backup.id]));
        };

        e_field = mat(metal.size(), 3);
        solve_efield_timing = timer::time([&]() {
            for (int i = 0; i < metal.size(); ++i) {
                const Id &id = metal[i];
                const auto [x_minus, x_plus, y_minus, y_plus, z_minus, z_plus] = geom.adjacent(id);

                e_field(i, span::all) = {
                        (V(x_minus, id) - V(x_plus, id)) / (2 * h),
                        (V(y_minus, id) - V(y_plus, id)) / (2 * h),
                        (V(z_minus, id) - V(z_plus, id)) / (2 * h)
                };
            }
        });

        field<std::string> header{"x", "y", "z", "Ex", "Ey", "Ez"};

        io_timing = timer::time([&]() {
            save(e_field, csv_name(join_string("../data/efield_", name, ".csv"), header));
            save_ovf(e_field, "../efield_data_text_rectangular.ovf");
        });
    }

    void solve_current() {
        const auto &[_, spins, metal, metal_idx] = geom;

        double rho_parallel = 10;
        double rho_perpendicular = 5;
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

        for (const auto [id, v]: applied_volts) {
            if (v == 0) {
                const double output = norm(current(*metal_idx[id.id], span::all));
                std::cout << id << " -> output = " << output << std::endl;
            }
        }

        if (save_files) {
            field<std::string> header{"x", "y", "z", "jx", "jy", "jz"};
            save(current, csv_name(join_string("../data/current_", name, ".csv"), header));
        }
    }

    void rotate_spin() {
        // perpendicular to s_hat so we can rotate around it
        vec3 rotation_axis;
        if (applied_volts.size() == 2) {
            rotation_axis = {s_hat(1), -s_hat(0), 0};
        } else if (applied_volts.size() == 3) {
            Id output = std::find_if(applied_volts.begin(), applied_volts.end(),
                                     [](const auto &pair) { return pair.second == 0; })->first;
            std::vector<vec3> vecs;
            for (const auto &[id, v]: applied_volts)
                if (v != 0)
                    vecs.emplace_back(conv_fixed<3, uword, double>(output.to_xyz() - id.to_xyz()));
            std::cout << "vecs[0] = " << normalise(vecs[0]) << std::endl;
            std::cout << "vecs[1] = " << normalise(vecs[1]) << std::endl;
            rotation_axis = cross(vecs[0], vecs[1]);
            rotation_axis = normalise(rotation_axis);
            double angle = acos(dot(vecs[0], vecs[1]) / (norm(vecs[0]) * norm(vecs[1])));
            std::cout << "angle = " << angle << std::endl;
            rotate_vector(rotation_axis, angle, s_hat).print("rotated");
        }
        rotation_axis = normalise(rotation_axis);
        std::cout << "rotation_axis = " << rotation_axis << std::endl;

        Id output_id = std::find_if(applied_volts.cbegin(), applied_volts.cend(),
                                    [](const auto &pair) { return pair.second == 0; })->first;
        uword output_idx = *geom.metal_idx[output_id.id];

        const int n_degrees = 360;
        const int steps_per_degree = 2;
        mat current_magnitude(n_degrees * steps_per_degree, 2);
        for (int degrees = 0; degrees < n_degrees * steps_per_degree; ++degrees) {
            const double deg = (double) degrees / steps_per_degree;
            const double theta = deg * datum::pi / 180.0;
            vec3 spin = rotate_vector(rotation_axis, theta, s_hat);
            std::fill(geom.spins.begin(), geom.spins.end(), spin);
            solve_current();
            rowvec out = current(output_idx, span::all);
            current_magnitude(degrees, span::all) = {deg, norm(out)};
        }
        field<std::string> header{"theta", "|j|"};
        current_magnitude.save(csv_name(join_string("../data/current_varied_", name, ".csv"), header));
    }
};

// 1 in 1 out
// vector n from in to out, rotate spin from negative n to positive n

int main() {
    /*namespace fs = std::filesystem;
    std::ofstream benchmark_output("../benchmark/swc1_laptop_superlu_csv.csv");
    benchmark_output << "size,sites,solve_linear_eq,solve_efield,io" << std::endl;

    for (const auto &entry: fs::directory_iterator("../benchmark")) {
        const auto &path = entry.path();
        if (path.extension() == ".ovf") {
            std::cout << path.stem().string() << std::endl;
            size_t points;

            Solver solver(path.stem().string(), "sw4c10.0_corners", SolverOpts{.use_closest_point = true});
            solver.solve_voltages();
            solver.solve_e_field();
            points = solver.geom.metal.size();

            benchmark_output << path.stem().string() << ","
                             << points << ","
                             << solver.solve_linear_eq_timing << ","
                             << solver.solve_efield_timing << ","
                             << solver.io_timing << std::endl;
        }
    }*/

//  Solver solver("1d", "1d",
//  Solver solver("2x2x2", "sw4c10.0_two_inputs",
//    for (const auto str: {"sw4c10.0_corners", "sw4c10.0_alt", "sw4c10.0_two_inputs"}) {
//        Solver solver("spin_wire4current_1_B=0.0", str,
//                      SolverOpts{.use_closest_point = true, .vary_spins = true});
//        solver.solve_voltages();
//        solver.solve_e_field();
////        solver.solve_current();
//        solver.rotate_spin();
//    }

    for (const auto str: {"sw4c10.0_two_inputs_alt", "sw4c10.0_more"}) {
//    for (const auto str: {"spin_wire4current_1_B=0.0", "spin_wire4current_2_B=0.0", "spin_wire4current_3_B=0.0"}) {
        Solver solver("spin_wire4current_1_B=0.0", str,
                      SolverOpts{.use_closest_point = true});
        solver.solve_voltages();
        solver.solve_e_field();
        solver.solve_current();
//        solver.rotate_spin();
    }

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
