#include "structures.h"
#include "obj_parser.h"

#include <unordered_set>

Flat getFlat(){
    Flat f;
    std::cin >> f.n.x >> f.n.y >> f.n.z;
    std::cin >> f.p.x >> f.p.y >> f.p.z;
    return f;
}

float PointInFlat (Vertex const &p, Flat const &f){//подставление координат точки в ур.плоскости
    return f.n.x * p.x + f.n.y * p.y + f.n.z * p.z + f.D();
}


void  PointClassify(Mesh &m , Flat const &f){//классификация точек тела относительно плоскости

    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < 0){//IN
            m.Vertices[i].c = 1;
        }
        else if( -0.0001 < PointInFlat(m.Vertices[i], f) < 0.0001){//ON
            m.Vertices[i].c = 0;
        }
        else if( PointInFlat(m.Vertices[i], f) > 0){//OUT
            m.Vertices[i].c = -1;
        }
    }
}

int InPoints(Mesh const &m , Flat const &f){//Количество точек in
    int in = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < 0){//IN
            in ++;
        }
    }
    return in;
}
int OnPoints(Mesh const &m , Flat const &f){//Количество точек on 
    int zero = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) == 0){//ON
            zero ++;
        }
    }
    return zero;
}
int OutPoints(Mesh const &m , Flat const &f){//Количество точек out
    int out = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) > 0){//OUT
            out ++;
        }
    }
    return out;
}

Vertex Segment_Flat_Intersection(Segment const &s, Flat const &f){//точка пересечения плоскости и отрезка
    if(PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0) {
        Vertex p = s.A - s.B; 
        
        float t = - ( PointInFlat(s.B, f)) / (PointInFlat(p, f) - f.D());

        p.x = (s.A.x * t) + (s.B.x * (1.f - t));
        p.y = (s.A.y * t) + (s.B.y * (1.f - t));
        p.z = (s.A.z * t) + (s.B.z * (1.f - t));
        return p;
    }
    throw std::runtime_error("Here should be intersection!!!");  //PointInFlat(p, f) - f.D() == 0 <=> A == B, но тогда невозможно, что PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0
}                                                                // эта функция вызывается для пар точек IN и OUT => PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0 всегда выполнено

int getVertexIndex(Vertex const &v, Mesh const &m){ //vector.find
    int n = m.Vertices.size();
    for(int i =0; i < m.Vertices.size(); ++i){
        n -= 1;
        if(operator==(m.Vertices[i], v ))
            return i;
    }
    return -1;
}

void DeleteVertex(Mesh &m, Vertex &v){
    std::vector<int> index_to_del;
    for(int i = 0; i < m.Vertices.size(); ++i)
        if( m.Vertices[i] == v){
            m.Vertices.erase(m.Vertices.begin() + i);
            return;
        }
            //index_to_del.push_back(i);
    //for(int j = 0; j < index_to_del.size(); ++j)
    //    m.Vertices.erase(m.Vertices.begin() + index_to_del[j] - j);
}

void DeleteFace(Mesh &m, Face &f){  //удаление грани из списка
    std::vector<int> index_to_del;
    for(int i = 0; i<m.Faces.size(); i++)
        if (m.Faces[i] == f)
            index_to_del.push_back(i);
    for(int j = 0; j < index_to_del.size(); ++j)
        m.Faces.erase(m.Faces.begin() + index_to_del[j] - j);
}

void DeleteIndexes(Mesh &m, Face &f, int code){
    std::vector<int> index_to_del;
    for(int i = 0; i < f.Indices.size(); ++i)
        if( m.Vertices[f.Indices[i]].c == code)
            index_to_del.push_back(i);
    for(int j = 0; j < index_to_del.size(); ++j)
        f.Indices.erase(f.Indices.begin() + index_to_del[j] - j); // проверить еще
}

void PushIndex(Face &f, int index, int j){//вставка индекса новой вершины в список на место j+1 (после j) в списке грани
    f.Indices.insert( std::next(f.Indices.begin(), j + 1) , index); 
}

bool isOnLine(Vertex &v, Vertex const &v1, Vertex const &v2){
    return Vector(v, v1).normalize() == Vector(v2, v).normalize();
}

std::vector<int> DuplicateVertecies(Mesh &m, std::vector<int> arr){//удаление по одинаковым индексам

    std::vector<int> d;
    std::vector<int> t;
    for(int i =0; i < arr.size() - 1; ++i)
        for(int j = i + 1; j < arr.size(); ++j){
            if(m.Vertices[arr[i]] == m.Vertices[arr[j]]){
                d.push_back(arr[j]);
                t.push_back(arr[i]);
            }     
        }

    for(int i = 0; i < d.size(); ++i){
        //DeleteVertex(m, m.Vertices[d[i]]);
        //std::cout << "deleted vertex:   "<< m.Vertices[d[i]].x <<' ' << m.Vertices[d[i]].y <<' ' << m.Vertices[d[i]].z<<' '<<std::endl;
        for(int n = 0; n < m.Faces.size(); ++n)
            for(int k = 0; k < m.Faces[n].Indices.size(); ++k){
                if(m.Faces[n].Indices[k] == d[i])
                    m.Faces[n].Indices[k] = t[i];
            }
    }

    return t;
}
    
void DeleteDuplicates(Mesh &m){
    std::vector<int> duplicates;
    for(int i = 0; i < m.Vertices.size() - 1; ++i)
        for(int j = i + 1; j < m.Vertices.size(); ++j)
            if(m.Vertices[i] == m.Vertices[j]){
                duplicates.push_back(j);
                std::cout << "dupl:  " << j << std::endl;
                std::cout << "ver  i:  " << i << ' '<<"u  "<< m.Vertices[i].x<<' '  <<  m.Vertices[i].y<<' ' <<  m.Vertices[i].z<<std::endl;
                std::cout << "ver  j:  " <<j <<' ' <<"u  "<< m.Vertices[j].x<<' '  <<  m.Vertices[j].y<<' ' <<  m.Vertices[j].z<<std::endl;
                std::cout << " divvvvvvvvvvv         "<< m.Vertices[i].x - m.Vertices[j].x <<' '<<   m.Vertices[i].y - m.Vertices[j].y <<' ' << m.Vertices[i].z - m.Vertices[j].z <<std::endl; 
            }
                
    for(int i = 0; i < duplicates.size(); ++i){
        std::cout << "sized   "<< duplicates[i] << std::endl;;
        DeleteVertex(m, m.Vertices[duplicates[i] - i]);
    }
}
            
    

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat const &f){//упорядочивание вершин в порядке обхода сечения
    std::vector<Vertex> tries;
    Vector norm = f.n.normalize();
    std::vector<int> del;

    std::unordered_set<int> unique_id_to_delete;
    for(int i = 0; i < intersect.size(); ++i)
        for(int j = 0; j < intersect.size()-1; ++j){
            if( j != i){
                for(int k = j+1; k < intersect.size(); ++k){
                    if(k != i && k != j){
                        if(isOnLine(intersect[i], intersect[j], intersect[k])){
                    //i между j и k
                            unique_id_to_delete.insert(i);
                            break;
                        }
                    }
                }
            }
        }
    
    std::vector<Vertex> new_vertices;
    for (int i = 0; i != intersect.size(); ++i) {
        if (!unique_id_to_delete.count(i))
            new_vertices.push_back(intersect[i]);
    }
    intersect = new_vertices;

    std::vector<int> arr(intersect.size()-1);//массив - количество появления вершин на второй позиции
    for(int i=0; i<intersect.size() - 1; ++i)
        arr[i] = 0;
    for( int i = 1; i < intersect.size() - 1; ++i ){
        for(int j = i + 1; j < intersect.size(); ++j ){

            Vector Try = Vector{intersect[0], intersect[i]}.cross(Vector{intersect[i], intersect[j]}).normalize();
            if (Try == norm)
                arr[i-1]++;
            else if (Try == -norm)
                arr[j-1]++;
        }
    }
    
        for (auto &&x: tries) {
        std::cout << x << '\n';
    }

    tries.push_back(intersect[0]);
    for(int j = 2; j <= intersect.size(); j++){ //tries - список вершин новой грани в порядке обхода. (сортировка)
        for(int i=0; i < intersect.size() - 1; ++i){
            if(arr[i] == intersect.size() - j)
                tries.push_back(intersect[i+1]);
        }
    }


    return tries;
}


void SpecialCases(Mesh &m, Flat const &f){
    int zero = OnPoints(m, f);//количество вершин объекта, лежащих на плоскости
    int sum = InPoints(m, f) - OutPoints(m, f);//cумма codes. нужно для определения крайних/вырожденных случаев

    if((sum == - (m.Vertices.size() - zero)) ){//случаи, где результат - пустое множество
        std::cerr << "Empty intersect " << std::endl;
        //можно отрисовать пустой экран
    }
    if((sum == (m.Vertices.size() - zero)) ){// случаи, когда результат - исходный объект
        std::cerr << "the object will not change " << std::endl;
        //отрисовка исходного меша
    }
}


Mesh ResultOfIntersect( Mesh const &m_in, Flat const &f){

    

    Mesh m{m_in};
    int index = m.Vertices.size();//индексы новых вершин
    Face new_face;
    std::vector<int> face_del;
    std::vector<Vertex> intersect; //points of intersect





    

    for(int i = 0; i < m.Faces.size(); ++i){//для каждой грани
        int out = 0;

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){//обходим каждую вершину грани

            if(m.Vertices[m.Faces[i].Indices[j]].c == -1){
                out++;
            }

            if( m.Vertices[m.Faces[i].Indices[j]].c  && m.Vertices[m.Faces[i].Indices[j]].c == -m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Segment s;
                s.A = m.Vertices[m.Faces[i].Indices[j]];
                s.B = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r = Segment_Flat_Intersection(s, f);
                r.c = 0;//ON
                m.Vertices.push_back(r);
                //new_face.Indices.push_back(index);
                
                PushIndex(m.Faces[i], index, j);// вставка значения index на место j
                index++;
            }
        }

        DeleteIndexes(m, m.Faces[i], -1);

        if(out == m.Faces[i].Indices.size()){
            face_del.push_back(i);
        }      

    }



    

    for(int i = 0; i < face_del.size(); ++i){
        DeleteFace(m, m.Faces[face_del[i] - i]);
    }
    int del = 0;
    for(int i = 0; i < m.Faces.size(); ++i){
        if(m.Faces[i].Indices.size() == 0){
            DeleteFace(m, m.Faces[i - del]);
            del++;
        }
    }


    

    std::vector<float> vert;

    //new_face.Indices = DuplicateVertecies(m, new_face.Indices);
    for(int i = 0; i < m.Vertices.size(); ++i){
        vert.push_back(m.Vertices[i].x);
        vert.push_back(m.Vertices[i].y);
        vert.push_back(m.Vertices[i].z);
    }



    DeleteDuplicates(m);

        
    for(int i =0; i < m.Vertices.size(); ++i){
        std::cout << "vert: "<< i <<"  ";
        std::cout <<  ' '<< m.Vertices[i].x << ' '<<  ' '<< m.Vertices[i].y <<  ' '<< m.Vertices[i].z << "  c:   "<<m.Vertices[i].c <<std:: endl;
    }
    std::cout << std:: endl;

    for(int i =0; i < m.Faces.size(); ++i){
        std::cout << "face: "<< i <<"  ";
        for(int j = 0; j < m.Faces[i].Indices.size(); ++j){
            std::cout <<  ' '<< m.Faces[i].Indices[j]<< ' ';
        }
        std::cout << std:: endl;
    }
    std::cout << std:: endl;



    std::vector<int> index_ver_del;//удаление вершин OUT
    for(int i = 0; i < m.Vertices.size(); ++i){
        if(m.Vertices[i].c == -1){
            index_ver_del.push_back(i);
        }
    }




    for(int j = 0; j < index_ver_del.size(); ++j){
        DeleteVertex(m, m.Vertices[index_ver_del[j] -j]);
    }



//////////////////////////////////////////////////////////////////
    for(int j =0; j < m.Faces.size(); ++j)
            
        for(int n =0; n < m.Faces[j].Indices.size(); ++n)
            for(int i = 0; i < vert.size(); i += 3){

                if(m.Faces[j].Indices[n] == i/3){
                    Vertex v = {{vert[i]}, {vert[i + 1]}, {vert[i + 2]}};
                    m.Faces[j].Indices[n] = getVertexIndex(v, m);
                    break;
                }
            }

    for(int i =0; i < m.Vertices.size(); ++i){//индексы сечения
        if(m.Vertices[i].c == 0){
            new_face.Indices.push_back(i);
        }
    }


    for (int i = 0; i < new_face.Indices.size(); ++i){
        intersect.push_back(m.Vertices[new_face.Indices[i]]);
    }

    intersect = tries(intersect, f);//вектор вершин в нужном порядке для новой грани

    Face face_intersect;
    for(int i =0; i < intersect.size(); ++i){
        face_intersect.Indices.push_back(getVertexIndex(intersect[i], m));
    }


    m.Faces.push_back(face_intersect);

    return m;
}

void Triangulation(Mesh &m){
    for(int i = 0; i < m.Faces.size(); ++i){
        if(m.Faces[i].Indices.size() > 3)
            
            for(int j = 0; j < m.Faces[i].Indices.size() - 1; j += 2){
                Face new_face;
                for( int n = 0; n < 3; ++n)
                    new_face.Indices.push_back(m.Faces[i].Indices[j + n]);
                m.Faces[i].Indices.erase(m.Faces[i].Indices.begin() + j + 1);
                m.Faces.push_back(new_face);
            }
    }
    /*
    for(int i =0; i < m.Vertices.size(); ++i){
        std::cout << "vert: "<< i <<"  ";
        std::cout <<  ' '<< m.Vertices[i].x << ' '<<  ' '<< m.Vertices[i].y <<  ' '<< m.Vertices[i].z <<std:: endl;
    }
    std::cout << std:: endl;

    for(int i =0; i < m.Faces.size(); ++i){
        std::cout << "face: "<< i <<"  ";
        for(int j = 0; j < m.Faces[i].Indices.size(); ++j){
            std::cout <<  ' '<< m.Faces[i].Indices[j]<< ' ';
        }
        std::cout << std:: endl;
    }
    std::cout << std:: endl;
    */

}