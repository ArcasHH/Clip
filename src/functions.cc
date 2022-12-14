#include "functions.h"

void Translation(Mesh &m, Vector const &r){
    for(int i =0; i<m.Vertices.size(); ++i){
        m.Vertices[i].x += r.x;
        m.Vertices[i].y += r.y;
        m.Vertices[i].z += r.z;
    }
}

void Scalation(Mesh &m, double x, double y, double z){
    for(int i =0; i<m.Vertices.size(); ++i){
        m.Vertices[i].x *= x;
        m.Vertices[i].y *= y;
        m.Vertices[i].z *= z;
    }
}

void Rotation(Mesh &m, double phi, double psi, double teta){
    
}