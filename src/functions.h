#include <iostream>
#include "structures.h"

void Translation(Mesh &m, Vector const &r);//сдвиг на вектор

void Scalation(Mesh &m, double x, double y, double z);//растяжение-сжатие по осям

void Rotation(Mesh &m, double phi, double psi, double teta);//вращения по углам Эйлера относительно 000