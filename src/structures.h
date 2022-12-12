#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>

#include "obj_parser.h"

struct Vector {
    double x{0.f}, y{0.f}, z{0.f};
    int c;

    Vector() = default;
    Vector(double x, double y, double z) : x{x}, y{y}, z{z} {}
    Vector(Vector const &A, Vector const &B) {
        x = B.x - A.x;
        y = B.y - A.y;
        z = B.z - A.z;
    }
    // Vector(Segment const &S) : Vector{S.A, S.B} {}

    double length_sq() const { return x*x + y*y + z*z; }
    double length() const { return std::sqrt(length_sq()); }
    Vector normalize() const { double const len = length(); return Vector{x / len, y / len, z / len};   }

    // -A; -> A.operator-()
    Vector operator-() {
        Vector tmp;
        tmp.x = -x;
        tmp.y = -y;
        tmp.z = -z;
        return tmp;
    }

    // B -= V -> B.operator(V)
    Vector &operator-=(Vector const &V) {
        x -= V.x;
        y -= V.y;
        z -= V.z;
        return *this;
    }

    Vector operator*(double t) const {
        Vector tmp;
        tmp.x *= t;
        tmp.y *= t;
        tmp.z *= t;
        return tmp;
    }

    Vector cross(Vector const &A) const{
        Vector r;
        r.x = (y * A.z) - (z * A.y);
        r.y = -(x * A.z) + (z * A.x);
        r.z = (x * A.y) - (y * A.x);
        return r;
    } //  it is VectorProduct. [this x A]
    double dot(Vector const &A) const{
        return x * A.x + y * A.y + z * A.z;
    } //  it is ScalarProduct. this * A
};

using Vertex = Vector;

inline Vector operator+(Vector const &B, Vector const &A) {
    Vector v;
    v.x = B.x + A.x;
    v.y = B.y + A.y;
    v.z = B.z + A.z;
    return v;
}

// C = A - B; -> A = .operator-(A, B)
inline Vector operator-(Vector const &B, Vector const &A) {
    Vector v{B};
    v -= A;
    return v;
}

inline bool operator==(Vector const &B, Vector const &A) {//с погрешностью в 1%
    double e = 1e-4;
    return (std::abs(A.x - B.x) < e && std::abs(A.y - B.y) < e && std::abs(A.z - B.z) < e);
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

inline bool operator==(Face const &f1, Face const &f2){
    if(f1.Indices.size() == f2.Indices.size()){
        for(int i = 0; i < f1.Indices.size(); ++i){
            for (int j = 0; j < f2.Indices.size(); ++j)
            if(f1.Indices[i] == f2.Indices[j])
                for(int k = 0; k < f1.Indices.size(); ++k)
                    if(f1.Indices[(i + k) % f1.Indices.size()] != f2.Indices[(j + k) % f2.Indices.size()])
                        return false;
        }
        return true;
    }
    return false;
}

struct Flat {
    Vector n;
    Vertex p;
    double D() const { return -n.x * p.x - n.y * p.y - n.z * p.z; }
};

inline bool operator==(Flat const &f1, Flat const &f2){
    if(f1.n == f2.n && f1.D() == f2.D())
        return true;
    else    
        return false;
}

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

double PointInFlat (Vertex const &p, Flat const &f);

int getVertexIndex(Vertex const &v, Mesh const &m);

void PointClassify(Mesh &m, Flat const &f);

void DeleteVertex(Mesh &m, Vertex &v);

void DeleteDuplicates(Mesh &m);

void DeleteIndexes(Mesh &m, Face &f, int code);

void DeleteFace(Mesh &m, Face &f);

void DeleteUncorrectFaces(Mesh &m);

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat const &f);

Vertex Segment_Flat_Intersection(Segment const &s, Flat const &f);

bool SpecialCases(Mesh &m, Flat const &f);

Mesh ResultOfIntersect( Mesh const &m_in, Flat const &f);

void Triangulation(Mesh &m);

bool Check(Mesh const &m);

void Intersect(Mesh &m, Flat const &f);

double cos(Vector const &v1, Vector const &v2);

