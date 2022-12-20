#pragma once

class Mesh;

Mesh getModel(char const * objFile);

void Write(Mesh const &m, double precise);