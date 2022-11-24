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

    while(std::getline(in, line)) {
        std::stringstream ss{line};
        char c;
        ss >> c;
        switch (c)
        {
        case 'v':
            Vertex v;
            ss >> v;
            pos.push_back(v);
            break;
        case 'f':
            Face f;
            ss >> f;
            faces.push_back(f);
            break;
        default:
            std::cerr << "unknown character: " << c << std::endl;
            break;
        }
    }

    return Mesh{pos, faces};
}