#ifndef MODEL_H
#define MODEL_H

#include <glad/glad.h> 

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <learnopengl/mesh.h>
#include <learnopengl/shader.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

Mesh loadOBJ(const char * path);
class Model 
{
public:
    // model data 
    vector<Texture> textures_loaded;	// stores all the textures loaded so far, optimization to make sure textures aren't loaded more than once.
    vector<Mesh>    meshes;
    string directory;
    bool gammaCorrection;

    // constructor, expects a filepath to a 3D model.
    Model(string const &path, bool gamma = false) : gammaCorrection(gamma)
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

Mesh loadOBJ(const char * path) {
	printf("Loading OBJ file %s...\n", path);

	std::vector<unsigned int> vertexIndices, uvIndices, normalIndices;
	std::vector<glm::vec3> temp_vertices; 
	std::vector<glm::vec2> temp_uvs;
	std::vector<glm::vec3> temp_normals;

    std::vector<Vertex> Vertices;


	FILE * file = fopen(path, "r");
	if( file == NULL ){
		throw std::runtime_error{"Impossible to open the file ! Are you in the right path ? See Tutorial 1 for details\n"};
	}

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
			std::string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
			int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2] );
			if (matches != 9){
				printf("File can't be read by our simple parser :-( Try exporting with other options\n");
				fclose(file);
				throw std::runtime_error{"File can't be read by our simple parser :-( Try exporting with other options\n"};
			}
			vertexIndices.push_back(vertexIndex[0]);
			vertexIndices.push_back(vertexIndex[1]);
			vertexIndices.push_back(vertexIndex[2]);
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

	// // For each vertex of each triangle
	// for( unsigned int i=0; i<vertexIndices.size(); i++ ){

	// 	// Get the indices of its attributes
	// 	unsigned int vertexIndex = vertexIndices[i];
	// 	unsigned int uvIndex = uvIndices[i];
	// 	unsigned int normalIndex = normalIndices[i];

    //     std::cout << "vI " << vertexIndex << '\n';
		
	// 	// Get the attributes thanks to the index
	// 	glm::vec3 vertex = temp_vertices[ vertexIndex-1 ];
	// 	glm::vec2 uv = temp_uvs[ uvIndex-1 ];
	// 	glm::vec3 normal = temp_normals[ normalIndex-1 ];
		
	// 	// Put the attributes in buffers
	// 	// out_vertices.push_back(vertex);
	// 	// out_uvs     .push_back(uv);
	// 	// out_normals .push_back(normal);

    //     // Vertices.push_back(Vertex{.Position = vertex, .Normal = normal, .TexCoords = uv});
	
	// }
	fclose(file);

    for (auto &&x : temp_vertices) {
        Vertices.push_back(Vertex{.Position = x});
    }

    for(auto &&x : vertexIndices) {
        --x;
    }
	
    return Mesh{Vertices, vertexIndices};
}


#endif
