#ifndef MODEL_H
#define MODEL_H

#include <glad/glad.h> 

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <learnopengl/mesh.h>
#include <learnopengl/shader_m.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
using namespace std;

namespace mdl{

Mesh loadOBJ(const char * path);
class Model 
{
public:
    vector<Mesh>    meshes;
    string directory;

    // constructor, expects a filepath to a 3D model.
    Model(string const &path)
    {
        loadModel(path);
    }

    // draws the model, and thus all its meshes
    void Draw(Shader &shader)
    {
        for(unsigned int i = 0; i < meshes.size(); i++)
            meshes[i].Draw(shader);
    }
    
private:
    // loads a model with supported ASSIMP extensions from file and stores the resulting meshes in the meshes vector.
    void loadModel(string const &path)
    {   
        meshes.push_back(loadOBJ(path.c_str()));
    }

};

glm::vec3 operator/(glm::vec3 const & v, int i) {
    glm::vec3 tmp = v;
    tmp.x /= i;
    tmp.y /= i;
    tmp.z /= i;
    return tmp;
}

Mesh loadOBJ(const char * path) {
	printf("Loading OBJ file %s...\n", path);

	std::vector<unsigned int> vertexIndices, uvIndices, normalIndices;
	std::vector<glm::vec3> temp_vertices; 
	std::vector<glm::vec2> temp_uvs;
	std::vector<glm::vec3> temp_normals;

    std::vector<glm::vec3> out_vertices;
	std::vector<glm::vec3> out_normals;

    std::vector<Vertex> Vertices;


	FILE * file = fopen(path, "r");
	if( file == NULL ){
		throw std::runtime_error{"Impossible to open the file ! Are you in the right path ? See Tutorial 1 for details\n"};
	}

    bool IsTriangulated = true;

	while( 1 ){

		char lineHeader[128];
		// read the first word of the line
		int res = fscanf(file, "%s", lineHeader);
		if (res == EOF)
			break; // EOF = End Of File. Quit the loop.

		// else : parse lineHeader
		
		if ( strcmp( lineHeader, "v" ) == 0 ){
			glm::vec3 vertex;
			fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
			temp_vertices.push_back(vertex);
		}else if ( strcmp( lineHeader, "vt" ) == 0 ){
			glm::vec2 uv;
			fscanf(file, "%f %f\n", &uv.x, &uv.y );
			uv.y = -uv.y; // Invert V coordinate since we will only use DDS texture, which are inverted. Remove if you want to use TGA or BMP loaders.
			temp_uvs.push_back(uv);
		}else if ( strcmp( lineHeader, "vn" ) == 0 ){
			glm::vec3 normal;
			fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
			temp_normals.push_back(normal);
		}else if ( strcmp( lineHeader, "f" ) == 0 ){
            unsigned int vertexIndex[4], uvIndex[3], normalIndex[3];
            int matches = fscanf(file, "%d %d %d %d", vertexIndex, vertexIndex+1, vertexIndex+2, vertexIndex+3);
            if (matches != 4) {
                matches = fscanf(file, "/%d/%d %d/%d/%d %d/%d/%d\n", &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2] );
                if (matches != 8){
                    printf("File can't be read by our simple parser :-( Try exporting with other options\n");
                    fclose(file);
                    throw std::runtime_error{"File can't be read by our simple parser :-( Try exporting with other options\n"};
                }
            } else {
                printf("faces\n");
                IsTriangulated = false;
            }
			vertexIndices.push_back(vertexIndex[0]);
			vertexIndices.push_back(vertexIndex[1]);
			vertexIndices.push_back(vertexIndex[2]);
            if (!IsTriangulated) vertexIndices.push_back(vertexIndex[3]);
			uvIndices    .push_back(uvIndex[0]);
			uvIndices    .push_back(uvIndex[1]);
			uvIndices    .push_back(uvIndex[2]);
			normalIndices.push_back(normalIndex[0]);
			normalIndices.push_back(normalIndex[1]);
			normalIndices.push_back(normalIndex[2]);
		}else{
			// Probably a comment, eat up the rest of the line
			char stupidBuffer[1000];
			fgets(stupidBuffer, 1000, file);
		}

	}

	fclose(file);

    glm::vec3 avgNV(temp_vertices.size());
    int count;

    std::vector<glm::vec3> vAndN(temp_vertices.size());

	for( unsigned int i=0; i<vertexIndices.size(); i++ ){

		// Get the indices of its attributes
		unsigned int vertexIndex = vertexIndices[i];
		unsigned int uvIndex = uvIndices[i];
		unsigned int normalIndex = normalIndices[i];
		
		// Get the attributes thanks to the index
		glm::vec3 vertex = temp_vertices[ vertexIndex-1 ];
		glm::vec2 uv = temp_uvs[ uvIndex-1 ];
		glm::vec3 normal = temp_normals[ normalIndex-1 ];
		
        vAndN[vertexIndex-1] += normal;
	}

    // for (auto &&x : temp_vertices) {
    //     Vertices.push_back(Vertex{.Position = x});
    // }

    for (size_t i = 0; i != temp_vertices.size(); ++i) {
        glm::vec3 N =  glm::normalize(vAndN[i]);
        Vertices.push_back(Vertex{.Position = temp_vertices[i], .Normal = N });
    }


    for(auto &&x : vertexIndices) {
        --x;
    }

    // std::vector<unsigned int> vi = {
    //     0,1,2,3,
    //     0,4,5,1,
    //     4,7,6,5,
    //     2,6,7,3,
    //     1,5,6,2,
    //     0,3,7,4
    // };
	
    return Mesh{Vertices, vertexIndices, IsTriangulated};
}

} // namespace mdl
#endif