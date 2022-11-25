#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>

#include "obj_parser.h"

struct Vector {
    float x{0.f}, y{0.f}, z{0.f};
    int c;

    Vector() = default;
    Vector(float x, float y, float z) : x{x}, y{y}, z{z} {}
    Vector(Vector const &A, Vector const &B) {
        x = B.x - A.x;
        y = B.y - A.y;
        z = B.z - A.z;
    }
    // Vector(Segment const &S) : Vector{S.A, S.B} {}

    float length_sq() const { return x*x + y*y + z*z; }
    float length() const { return std::sqrt(length_sq()); }
    Vector normalize() const { float const len = length(); return Vector{x / len, y / len, z / len};   }

    // -A; -> A.operator-()
    Vector &operator-() {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }

    // B -= V -> B.operator(V)
    Vector &operator-=(Vector const &V) {
        x -= V.x;
        y -= V.y;
        z -= V.z;
        return *this;
    }

    Vector operator*(float t) const {
        Vector tmp;
        tmp.x *= t;
        tmp.y *= t;
        tmp.z *= t;
        return tmp;
    }

    Vector cross(Vector const &A) const; // TODO it is VectorProduct. [this x A]
    float dot(Vector const &A) const; // TODO it is ScalarProduct. this * A
};
using Vertex = Vector;

inline Vector operator+(Vector const &B, Vector const &A) {
    Vector v;
    v.x = B.x + A.x;
    v.y = B.y + A.y;
    v.z = B.z + A.z;
    return v;
}
inline bool operator==(Vector const &B, Vector const &A) {
    return B.x == A.x && B.y == A.y && B.z == A.z;
}

// C = A - B; -> A = .operator-(A, B)
inline Vector operator-(Vector const &B, Vector const &A) {
    Vector v{B};
    v -= A;
    return v;
}

struct Segment {
    Vertex A;
    Vertex B;
    Segment() = default;
    Segment(Vertex const &A, Vertex const &B) : A{A}, B{B} {}
};


struct Face {
    std::vector<int> Indices;////////////////////////////пересмотреть
    
    int operator[] (int pos) const { return Indices.at(pos); }
    int &operator[] (int pos) { return Indices.at(pos); }
};

struct Flat {
    Vector n;
    Vertex p;
    float &A = n.x;
    float &B = n.y;
    float &C = n.z;

    // TODO: Do or A(), B() ... or n.x, n.y...
    float D() const { return -A * p.x - B * p.y - C * p.z; }
};

inline std::istream &operator>>(std::istream &is, Face &f) {
    is >> f[0] >> f[1] >> f[2] >> f[3];
    return is;
}
inline std::istream &operator>>(std::istream &is, Vertex &v) {
    is >> v.x >> v.y >> v.z;
    return is;
}
inline std::ostream &operator<<(std::ostream &os, Face &f) {
    os << f[0] << ' ' <<  f[1] << ' ' <<  f[2] << ' ' <<  f[3];
    return os;
}
inline std::ostream &operator<<(std::ostream &os, Vertex &v) {
    os << v.x << ' ' << v.y << ' ' << v.z;
    return os;
}


// Vertex representation
struct Mesh {
    std::vector<Vertex> Vertices;
    std::vector<Face> Faces;

    void dump() const {
        for (Vertex x : Vertices) {
            std::cout << "v " << x << std::endl;
        }
        for (Face x : Faces) {
            std::cout << "f " << x << std::endl;
        }
    }
};

float ScalarProduct(Vector &v1, Vector &v2);

Vector VectorProduct(Vector &v1, Vector &v2);

int PointInFlat (Vertex const &p, Flat const &f);

void PointClassify(Vertex &p, Flat &f);

void DeleteVertex(Mesh &m, Vertex &v);

void DeleteIndexes(Mesh &m, Face &f, int code);

void PushIndex(Face &f, int index);

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat f);

Vertex Segment_Flat_Intersection(Segment const &s, Flat const &f);

Mesh ResultOfIntersect( Mesh &m, Flat &f);
