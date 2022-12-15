#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "structures.h"

Mesh getModel(char const * objFile)
{
    std::ifstream in{objFile};
    if(!in.is_open())
        throw std::invalid_argument{"can not find file " + std::string{objFile}};

    std::vector<Vertex> pos;
    std::vector<Face> faces;

    std::string line; 
    ///jg



    return Mesh{pos, faces};
}