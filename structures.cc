#include "structures.h"
#include "obj_parser.h"


float Vector::dot(Vector const &A) const {
    return x * A.x + y * A.y + z * A.z;
}

Vector Vector::cross(Vector const &A) const {
    Vector r;
    r.x = (y * A.z) - (z * A.y);
    r.y = -(x * A.z) + (z * A.x);
    r.z = (x * A.y) - (y * A.x);
    return r;
}

int PointInFlat (Vertex const &p, Flat const &f){//подставление координат точки в ур.плоскости
    return f.A * p.x + f.B * p.y + f.C * p.z + f.D();
}


void  PointClassify(Mesh &m , Flat &f){//классификация точек тела относительно плоскости

    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < 0){//IN
            m.Vertices[i].c = 1;
        }
        else if( PointInFlat(m.Vertices[i], f) == 0){//ON
            m.Vertices[i].c = 0;
        }
        else if( PointInFlat(m.Vertices[i], f) > 0){//OUT
            m.Vertices[i].c = -1;
        }
    }
}

int InPoints(Mesh &m , Flat &f){//Количество точек in
    int in = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < 0){//IN
            in ++;
        }
    }
    return in;
}
int OnPoints(Mesh &m , Flat &f){//Количество точек on 
    int zero = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) == 0){//ON
            zero ++;
        }
    }
    return zero;
}
int OutPoints(Mesh &m , Flat &f){//Количество точек out
    int out = 0;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) > 0){//OUT
            out ++;
        }
    }
    return out;
}

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat f){//упорядочивание вершин в порядке обхода сечения
    std::vector<Vertex> tries;
    Vector norm = f.n.normalize();

    int num_of_triangles = 0; //количество треугольников, на которые разбивается сечение
    for(int i = intersect.size() - 2; i>0; i--){
        num_of_triangles +=i;
    }

    std::vector<int> arr(intersect.size()-1);//массив - количество появления вершин на второй позиции
    for(int i=0; i<intersect.size() - 1; ++i)
        arr[i] = 0;
    for( int i = 1; i < intersect.size() - 2; ++i ){
        for(int j = i + 1; j < intersect.size(); ++j ){

            Vector Try = Vector{intersect[0], intersect[i]}.cross(Vector{intersect[i], intersect[j]}).normalize();
            
            if (Try == norm)
                arr[i-1]++;
            else if (Try == -norm)
                arr[j-1]++;
        }
    }

    for(int j = 2; j <= intersect.size(); j++){ //tries - список вершин новой грани в порядке обхода. (сортировка)
        for(int i=0; i < intersect.size() - 1; ++i){
            if(arr[i] == intersect.size() - j)
                tries.push_back(intersect[i+1]);
        }
    }
    return tries;
}

Vertex Segment_Flat_Intersection(Segment const &s, Flat const &f){//точка пересечения плоскости и отрезка
    if(PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0) {
        Vertex p = s.A - s.B;
        float t = - ( PointInFlat(s.B, f)) / (PointInFlat(p, f) - f.D());
        return (s.A * t) + (s.B * (1.f - t));
    }

    throw std::runtime_error("Here should be intersection!!!");
}

int getVertexIndex(Vertex v, Mesh m){
    int n = m.Vertices.size();
    for(int i =0; i < m.Vertices.size(); ++i){
        n = n-1;
        if(operator==(m.Vertices[i], v ))
            return i;
    }
    return -1;
}

void DeleteVertex(Mesh &m, Vertex &v){
    for(int i =0; i < m.Vertices.size(); ++i){ //удаление 1 вершины из списка объекта//////////////////////////////
        if( m.Vertices[i] == v){//////////////////////////////////////
            if (i == m.Vertices.size() - 1){
                m.Vertices.pop_back();
                return;
            }
            else {
                for(int j = i + 1; j < m.Vertices.size() ; ++j){
                    m.Vertices[j - 1] = m.Vertices[ j ];
                }
            }
        }
    }
}

void DeleteIndexes(Mesh &m, Face &f, int code){
    for(int i =0; i < f.Indices.size(); ++i){ //удаление из списка вершины//////////////////////////////
            if( m.Vertices[f.Indices[i]].c == code){
                for(int j = i + 1; j < f.Indices.size() ; ++j){
                    f.Indices[j - 1] = f.Indices[j];
                }
                f.Indices.pop_back();
            }
        }
}

void PushIndex(Face &f, int index, int j){
    //вставка индекса index новой вершины в конец списвка грани
    if(j == f.Indices.size() - 1){
        f.Indices.push_back(index);
    }
    else{
    //вставка индекса новой вершины в список на место j+1 (после j) в списке грани
        f.Indices.push_back( f.Indices[f.Indices.size() - 1] );//копируем последний элемент в конец
        for(int i = f.Indices.size() - 2; i > j; i--){
            f.Indices[i + 1] = f.Indices[i];
        }
        f.Indices[j+1] = index;
    }
}


//не работает, если плоскость пересекает вершины корректно

Mesh ResultOfIntersect( Mesh const &m_in, Flat &f){
    Mesh m{m_in};

    int index = m.Vertices.size();//индексы новых вершин
    int zero = OnPoints(m, f);//количество вершин объекта, лежащих на плоскости
    int sum = InPoints(m, f) - OutPoints(m, f);//cумма codes. нужно для определения крайних/вырожденных случаев

    Face new_face;
    m.Faces.push_back(new_face);
    
    if((sum == - (m.Vertices.size() - zero)) ){//случаи, где результат - пустое множество
        std::cerr << "Empty intersect " << std::endl;
        return m;
    }
    if((sum == (m.Vertices.size() - zero)) ){// случаи, когда результат - исходный объект
        return m;
    }

    std::vector<Vertex> intersect; //points of intersect
    for(int i = 0; i < m.Faces.size(); ++i){//для каждой грани
        int out = 0;

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){//обходим каждую вершину грани
            if(m.Vertices[m.Faces[i].Indices[j]].c == -1){
                out++;
            }

            if( m.Vertices[m.Faces[i].Indices[j]].c && m.Vertices[m.Faces[i].Indices[j]].c == -m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Segment s;
                s.A = m.Vertices[m.Faces[i].Indices[j]];
                s.B = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r = Segment_Flat_Intersection(s, f);
                r.c = 0;//ON
                intersect.push_back(r);
                m.Vertices.push_back(r);
                
                PushIndex(m.Faces[i], index, j);
                index++;
            }
        }

        DeleteIndexes(m, m.Faces[i], -1);

        

    }
    // for(int i = 0; i < m.Faces.size(); ++i){
    //     for(int j =0; j < m.Faces[i].Indices.size(); ++j)
    //         if(out == 4){//удаление i-ой грани из списка
    //             for(int j = i+1; j<m.Faces.size(); j++){
    //                 m.Faces[j-1] = m.Faces[j];
    //         }
    // }
        
    
    // intersect = tries(intersect, f);//вектор вершин в нужном порядке для новой грани
    // for(int i = 0; i<intersect.size(); ++i){
    //     new_face.Indices.push_back(getVertexIndex(intersect[i], m));
    // }
    
    // m.Faces.push_back(new_face);

    // }
    return m;
}