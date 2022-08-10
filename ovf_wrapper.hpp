#ifndef MODELING_OVF_WRAPPER_HPP
#define MODELING_OVF_WRAPPER_HPP

#include <array>
#include <vector>
#include <memory>
#include <armadillo>
#include "libraries/ovf/include/ovf.h"

namespace ovf {

struct MeshType {
    enum Value {
        Rectangular,
        Irregular
    };

    MeshType() = delete;

    constexpr MeshType(Value v) : value(v) {}

    constexpr operator Value() const {
        return value;
    }

    explicit constexpr operator char *() const {
        switch (this->value) {
            case Rectangular:
                return const_cast<char *>("rectangular");
            case Irregular:
                return const_cast<char *>("irregular");
            default:
                throw std::runtime_error("Unknown mesh type");
        };
    }

    explicit operator bool() const = delete;

private:
    Value value;
};

struct Segment {
    std::unique_ptr<ovf_segment> ptr;

    ovf_segment *operator->();

    const ovf_segment *operator->() const;

    Segment() : ptr(std::unique_ptr<ovf_segment>(ovf_segment_create())) {}

    Segment title(const char *title) {
        ptr->title = const_cast<char *>(title);
        return std::move(*this);
    }

    Segment description(const char *description) {
        ptr->comment = const_cast<char *>(description);
        return std::move(*this);
    }

    Segment value_dim(int value_dim) {
        ptr->valuedim = value_dim;
        return std::move(*this);
    }

    Segment value_units(const char *value_units) {
        ptr->valueunits = const_cast<char *>(value_units);
        return std::move(*this);
    }

    Segment value_labels(const char *value_labels) {
        ptr->valuelabels = const_cast<char *>(value_labels);
        return std::move(*this);
    }

    Segment mesh_unit(const char *mesh_unit) {
        ptr->meshunit = const_cast<char *>(mesh_unit);
        return std::move(*this);
    }

    Segment bounds(std::array<float, 3> bounds_min, std::array<float, 3> bounds_max) {
        for (int i = 0; i < 3; ++i) {
            ptr->bounds_min[i] = bounds_min[i];
            ptr->bounds_max[i] = bounds_max[i];
        }
        return std::move(*this);
    }

    Segment mesh_type(MeshType mesh) {
        ptr->meshtype = static_cast<char *>(mesh);
        return std::move(*this);
    }

    Segment rectangular_mesh() {
        return mesh_type(MeshType::Rectangular);
    }

    Segment irregular_mesh() {
        return mesh_type(MeshType::Irregular);
    }

    bool is_rectangular() const {
        return strcmp(ptr->meshtype, "rectangular") == 0;
    }

    bool is_irregular() const {
        return strcmp(ptr->meshtype, "irregular") == 0;
    }

    /// only valid on rectangular meshes
    Segment n_nodes(std::array<int, 3> n_nodes) {
        int N = 1;
        for (int i = 0; i < 3; ++i) {
            int n = n_nodes[i];
            ptr->n_cells[i] = n;
            N *= n;
        }
        ptr->N = N;

        return std::move(*this);
    }

    /// only valid on rectangular meshes
    Segment step_size(std::array<float, 3> step_size) {
        for (int i = 0; i < 3; ++i) {
            ptr->step_size[i] = step_size[i];
        }
        return std::move(*this);
    }

    /// only valid on rectangular meshes
    Segment base(std::array<float, 3> origin) {
        for (int i = 0; i < 3; ++i) {
            ptr->origin[i] = origin[i];
        }
        return std::move(*this);
    }

    /// only valid on irregular meshes
    Segment point_count(int point_count) {
        ptr->pointcount = point_count;
        return std::move(*this);
    }

    /// Makes a vector with the exact size to hold all the mesh information in this segment.
    template<typename T>
    std::vector<T> allocate_data_vector(
            typename std::enable_if<std::is_floating_point<T>::value>::type * = 0
    ) {
        if (is_rectangular()) {
            return std::vector<T>(ptr->valuedim * ptr->N);
        } else if (is_irregular()) {
            return std::vector<T>((ptr->valuedim + 3) * ptr->pointcount);
        } else {
            throw std::runtime_error("unreachable");
        }
    }

    /// Makes a matrix with the exact size to hold all mesh information in this segment.
    ///
    /// Note: due to how arma matrices are stored in memory, this matrix is a transpose of the "actual" shape you want
    template<typename T>
    arma::Mat<T> allocate_data_mat(
            typename std::enable_if_t<std::is_floating_point_v<T>> * = 0
    ) {
        if (is_rectangular()) {
            return arma::Mat<T>(ptr->valuedim, ptr->N);
        } else if (is_irregular()) {
            return arma::Mat<T>(ptr->valuedim + 3, ptr->pointcount);
        } else {
            throw std::runtime_error("unreachable");
        }
    }
};

enum Format {
    Text,
    Binary
};

struct File {
    ovf_file *ptr;
    std::vector<Segment> segments;
    /// Whether to throw an error or just print the error message. Defaults to crashing.
    bool fatal_errors = true;

    ovf_file *operator->() const;

    explicit File(ovf_file *ptr) : ptr(ptr) {}

    static File open_file(const char *path);

    static File create_file(const char *path);

    ~File();

    [[nodiscard]] arma::mat read(int idx = 0);

    [[nodiscard]] std::vector<arma::mat> read_all();

    [[nodiscard]] static arma::mat read_file(const char *path, int idx = 0);

    void write(Segment segment, const double *data, Format format = Format::Binary);

    static void write_file(const char *path, Segment segment, const double *data, Format format = Format::Binary);

private:
    template<typename F, typename ...Args>
    void handle_err(const char *ctx, F ovf_function, Args... args) {
        int err = ovf_function(ptr, args...);
        if (err != OVF_OK) {
            std::ostringstream msg;
            if (ctx) {
                msg << ctx << ": ";
            }
            msg << ovf_latest_message(ptr);
            if (fatal_errors) {
                throw std::runtime_error(msg.str());
            } else {
                std::cout << msg.str() << std::endl;
            }
        }
    };
};

}
#endif //MODELING_OVF_WRAPPER_HPP
