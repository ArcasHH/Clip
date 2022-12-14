#include <iostream>
#include "structures.h"

Flat FlatByPoints(Vertex const&v1, Vertex const&v2, Vertex const&v3);

void Test(Mesh const &m_in, double precise);

void Icosahedron( Mesh &m, double precise);
void Cuboctahedron(Mesh &m, double precise);

void Rhombicuboctahedron(Mesh &m, double precise);
void Rhombicuboctahedron2(Mesh &m, double precise);
void Rhombicuboctahedron3(Mesh &m, double precise);

void Pyramid(Mesh &m, double precise);

void Octahedron(Mesh &m, double precise);

void Tetrahedron(Mesh &m, double precise);