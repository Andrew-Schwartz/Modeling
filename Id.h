#ifndef MODELING_ID_H
#define MODELING_ID_H

#include <armadillo>

using namespace arma;

struct Id {
    static uword SizeX;
    static uword SizeY;
    static uword SizeZ;

    static std::array<uword, 3> Sizes() {
        return {SizeX, SizeY, SizeZ};
    }

    uword id;

    explicit Id(uword id);

    Id(uword x, uword y, uword z);

    explicit Id(const uvec3 &xyz);

    Id operator+(uword other) const;

    Id operator+(Id other) const;

    /// returns std::nullopt if this operation would result in an Id outside of the sample
    std::optional<Id> operator+(const ivec3 &displacement) const;

    Id operator-(uword other) const;

    Id operator-(Id other) const;

    bool operator==(const Id &rhs) const;

    bool operator!=(const Id &rhs) const;

    bool operator<(const Id &rhs) const;

    bool operator>(const Id &rhs) const;

    bool operator<=(const Id &rhs) const;

    bool operator>=(const Id &rhs) const;

    friend std::ostream &operator<<(std::ostream &os, const Id &id);

    [[nodiscard]] uvec3 to_xyz() const;

    [[nodiscard]] double dist3d(const Id &other) const;
};

#endif //MODELING_ID_H
