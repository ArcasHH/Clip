#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "structures.h"

void Write(Mesh const &m){
    std::ofstream fout("C:\\Users\\Arcasha\\Desktop\\tmp\\cut\\Clip\\models\\new.obj");
    //fout.open("cppstudio.txt");
    fout << "# Blender3D v249 OBJ File: untitled.blend" <<std::endl;
    fout<<"# www.blender3d.org" <<std::endl;
    fout<<"mtllib cube.mtl" <<std::endl;

    for(int i = 0; i < m.Vertices.size(); ++i){
        fout << "v" <<' '<< m.Vertices[i].x <<' '<< m.Vertices[i].y <<' '<< m.Vertices[i].z <<' '<< std::endl;
    }

    fout<<"vt"<<' '<< 0.748573 << ' '<< 0.750412 << std::endl;

    std::vector<Vector> normals;
    for(int i =0; i < m.Faces.size(); ++i){
        Vector n = Vector{{m.Vertices[m.Faces[i].Indices[0]]},{m.Vertices[m.Faces[i].Indices[1]]}}.cross(Vector{{m.Vertices[m.Faces[i].Indices[1]]},{m.Vertices[m.Faces[i].Indices[2]]}});
        bool uniq = true;
        if(n.length() > 0){
            n = n.normalize();
            if(normals.size() != 0)
                for(int j =0; j < normals.size(); ++j)
                    if(normals[j] == n)
                        uniq = false; 
            if(uniq){
                normals.push_back(n);
                fout<<"vn"<< ' '<< n.x <<' '<< n.y<<' '<< n.z<<std::endl;
            }   
        }
    }

    fout<<"usemtl Material_ray.png"<<std::endl;
    fout<<"s off"<<std::endl;
    int normal_index = 1;
    for(int i =0; i < m.Faces.size(); ++i){
        Vector n = Vector{{m.Vertices[m.Faces[i].Indices[0]]},{m.Vertices[m.Faces[i].Indices[1]]}}.cross(Vector{{m.Vertices[m.Faces[i].Indices[1]]},{m.Vertices[m.Faces[i].Indices[2]]}}); //normal
        //std::cout<< n <<std::endl;
        for(int j =0; j < normals.size(); ++j){
            if( n.length() > 0){
                n = n.normalize();
                if( n == normals[j])
                    normal_index = j + 1;
            }
                
        }
        fout<<"f"<<' ';
        for(int j =0; j < m.Faces[i].Indices.size(); ++j){
            fout<< m.Faces[i].Indices[j] + 1<<"/1/"<<normal_index<<' ';
        }
        fout<<std::endl;
    }




    //system("pause");
    fout.close();
    return;
}


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
        {
            Vertex v;
            ss >> v;
            pos.push_back(v);
            break;
        }
        case 'f':
        {
            Face f;
            ss >> f;
            faces.push_back(f);
            break;
        }
        default:
            std::cerr << "unknown character: " << c << std::endl;
            break;
        }
    }

    return Mesh{pos, faces};
}