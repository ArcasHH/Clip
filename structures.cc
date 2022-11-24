#include "structures.h"
#include "obj_parser.h"


float ScalarProduct(Vector &r1, Vector &r2){
    return (r1.x * r2.x + r1.y * r2.y + r1.z * r2.z);
};

struct Vector VectorProduct(Vector r1, Vector r2){
    struct Vector r;
    r.x = (r1.y * r2.z) - (r1.z * r2.y);
    r.y = -(r1.x * r2.z) + (r1.z * r2.x);
    r.z = (r1.x * r2.y) - (r1.y * r2.x);
    return r;
};

int PointInFlat (Vertex &p, Flat &f){
    return f.A * p.x + f.B * p.y + f.C * p.z + f.D;
}

int  PointClassify(Mesh &m , Flat &f){
    int zero = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) > 0){ //IN
            m.Vertices[i].c = 1;
        }
        else if( PointInFlat(m.Vertices[i], f) == 0){//ON
            m.Vertices[i].c = 0;
            zero ++;
        }
        else if( PointInFlat(m.Vertices[i], f) < 0){//OUT
            m.Vertices[i].c = -1;
        }
    }
    return zero;
}

Face NewFace(std::vector<Vertex> &intersect, Flat f){
    Face new_face;
 // построение выпуклой оболочки
    if(intersect.size() == 3){
        struct Vector norm = f.n.normalize();
        struct Vector try = VectorProduct( Vector::Vector(intersect[0], intersect[1]), Vector::Vector(intersect[1], intersect[2]));
        if(try.x == norm.x && try.y == norm.y && try.z == norm.z ){
            new_face.Indices.push_back(2, 1, 0);
        }
        else if(try.x == -norm.x && try.y == -norm.y && try.z == -norm.z )
            new_face.Indices.push_back(0, 1, 2);
    }
}

Mesh ResultOfIntersect( Mesh &m, Flat &f){
    int index = m.Vertices.size();
    int zero = PointClassify(m, f);
    int sum;
    for(int i =0; i<m.Vertices.size(); ++i){
        sum += m.Vertices[i].c;
    }
    if(sum == - (m.Vertices.size() - zero)){
        std::cerr << "Empty intersect " << std::endl;
        return;
    }
    if(sum == (m.Vertices.size() - zero)){
        return m;
    }

    std::vector<Vertex> intersect; //points of intersect
    for(int i = 0; i < m.Faces.size(); ++i){

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){

            if( m.Vertices[m.Faces[i].Indices[j]].c != m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Vertex p = m.Vertices[m.Faces[i].Indices[j]];
                Vertex q = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r;
                float t = - (f.A * q.x + f.B * q.y + f.C * q.z + f.D) / (f.A * (p.x - q.x) + f.B * (p.y - q.y) + f.C * (p.z - q.z));
                r.x = p.x * t + q.x * (1 - t);
                r.y = p.y * t + q.y * (1 - t);                       //можно упрость написав операторы умножения на число и сложения для точек
                r.z = p.z * t + q.z * (1 - t);
                r.c = 0;
                intersect.push_back(r);
                index ++;

                //вставка вершины в конец списвка грани
                if(j == m.Faces[i].Indices.size() - 1){
                    m.Faces[i].Indices.push_back(index);
                }
                else{
                //вставка вершины в список на место j+1 (после j) в списке грани
                    m.Faces[i].Indices.push_back( m.Faces[i].Indices[m.Faces[i].Indices.size() - 1] );//копируем последний элемент в конец
                    for(int s = m.Faces[i].Indices.size() - 2; s > j; s--){
                        m.Faces[i].Indices[s] = m.Faces[i].Indices[s-1];
                    }
                    m.Faces[i].Indices[j+1] = index;
                }
            }

        }

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){ //удаление из списка вершины//////////////////////////////
            if( m.Vertices[m.Faces[i].Indices[j]].c == -1){
                for(int s = j; s < m.Faces[i].Indices.size() - 1 ; ++s){
                    m.Faces[i].Indices[s] = m.Faces[i].Indices[s+1];
                }
                m.Faces[i].Indices.pop_back;
            }
        }

    }
    Face intersect_face = NewFace(intersect, f);


}