#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <iostream>

#include "obj_parser.h"

struct Vertex {
    float x, y, z;
    int c;
};

struct Segment {
    Vertex A;
    Vertex B;
    Segment(Vertex A, Vertex B) : A{A}, B{B} {}
};

struct Vector {
    float x, y, z;

    Vector(float x, float y, float z) : x{x}, y{y}, z{z} {}
    Vector(Vertex A, Vertex B) {
        x = B.x - A.x;
        y = B.y - A.y;
        z = B.z - A.z;
    }
    Vector(Segment const &S) : Vector{S.A, S.B} {}

    float length_sq() const { return x*x + y*y + z*z; }
    float length() const { return std::sqrt(length_sq()); }
    Vector normalize() const { float const len = length(); return Vector{x / len, y / len, z / len};   }
};

struct Face {
    std::vector<int> Indices;////////////////////////////пересмотреть
    
    int operator[] (int pos) const { return Indices.at(pos); }
    int &operator[] (int pos) { return Indices.at(pos); }
};

struct Flat {
    Vector n;
    Vertex p;
    float A = n.x;
    float B = n.y;
    float C = n.z;
    float D = -A * p.x - B * p.y - C * p.z;
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

Face NewFace(std::vector<Vertex> &intersect, Flat f);

int PointInFlat (Vertex &p, Flat &f);

void PointClassify(Vertex &p, Flat &f);

Mesh ResultOfIntersect( Mesh &m, Flat &f);

