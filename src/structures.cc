#include "structures.h"
#include "obj_parser.h"

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
        else if( PointInFlat(m.Vertices[i], f) == 0){//ON
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
        return (s.A * t) + (s.B * (1.f - t));
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
    for(int i =0; i < m.Vertices.size(); ++i) //удаление 1 вершины из списка объекта//////////////////////////////
        if( m.Vertices[i] == v){
             m.Vertices.erase(m.Vertices.begin() + i - 1);
            /*
            if (i == m.Vertices.size() - 1){
                m.Vertices.pop_back();
                return;
            }
            else 
                for(int j = i + 1; j < m.Vertices.size() ; ++j)
                    m.Vertices[j - 1] = m.Vertices[ j ];
            */

        }
}

void DeleteFace(Mesh &m, Face &f){  //удаление грани из списка
for(int i = 0; i<m.Faces.size(); i++)
    if (m.Faces[i] == f){
        m.Faces.erase(m.Faces.begin() + i - 1);
        /*
        if (i == m.Faces.size() - 1){
            m.Faces.pop_back();
            return;
        }
        else 
            for(int j = i + 1; j < m.Faces.size() ; ++j)
                m.Faces[j - 1] = m.Faces[ j ];
        */
    }
}

void DeleteIndexes(Mesh &m, Face &f, int code){
    for(int i =0; i < f.Indices.size(); ++i) //удаление вершины из списка грани
        if( m.Vertices[f.Indices[i]].c == code)
                f.Indices.erase(f.Indices.begin() + i - 1);
}


void PushIndex(Face &f, int index, int j){//вставка индекса новой вершины в список на место j+1 (после j) в списке грани
    f.Indices.insert( std::next(f.Indices.begin(), j + 1) , index);
/*
//вставка индекса index новой вершины в конец списвка грани
    if(j == f.Indices.size() - 1){
        f.Indices.push_back(index);
    }
    else{
    
        f.Indices.push_back( f.Indices[f.Indices.size() - 1] );//копируем последний элемент в конец
        for(int i = f.Indices.size() - 2; i > j; i--){
            f.Indices[i + 1] = f.Indices[i];
        }
        f.Indices[j+1] = index;
    }
*/
    
}

bool isOnLine(Vertex &v, Vertex const &v1, Vertex const &v2){
    if(Vector(v, v1).normalize() == Vector(v2, v1).normalize())
        return true;
    return false;
}

void DeleteDuplicateVertecies(Mesh &m, std::vector<int> arr){
    std::vector<int> d;
    for(int i =0; i < arr.size() - 1; ++i)
        for(int j = i + 1; j < arr.size(); ++j)

            if(m.Vertices[arr[i]] == m.Vertices[arr[j]]){ //сравнение с погрешностью
                DeleteVertex(m, m.Vertices[arr[j]]);
                d.push_back(arr[j]);
                for(int j =0; j < m.Faces.size(); ++j)
                    for( int n = 0; n <3; ++n)
                        if(m.Faces[j].Indices[n] == arr[j])
                            m.Faces[j].Indices[n] = arr[i];
            }
    
    for(int i =0; i < d.size() - 1; ++i)
        for(int j = 0; j < arr.size(); ++j)
            if (d[i] == arr[j])
                arr.erase(arr.begin() + j - 1);
}

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat const &f){//упорядочивание вершин в порядке обхода сечения
    std::vector<Vertex> tries;
    Vector norm = f.n.normalize();

    for(int i = 1; i < intersect.size() + 1; ++i){
        if(isOnLine(intersect[i % intersect.size()], intersect[i-1], intersect[(i+1) % intersect.size()])){
            intersect.erase(intersect.begin() + i - 1);
        }
    }

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

    std::vector<Vertex> intersect; //points of intersect
    for(int i = 0; i < m.Faces.size(); ++i){//для каждой грани
        int out = 0;

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){//обходим каждую вершину грани

            if(m.Vertices[m.Faces[i].Indices[j]].c == -1){
                out++;
            }
            

            if(m.Vertices[m.Faces[i].Indices[j]].c == 0){
                intersect.push_back(m.Vertices[m.Faces[i].Indices[j]]);
            }

            else if( m.Vertices[m.Faces[i].Indices[j]].c  && m.Vertices[m.Faces[i].Indices[j]].c == -m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Segment s;
                s.A = m.Vertices[m.Faces[i].Indices[j]];
                s.B = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r = Segment_Flat_Intersection(s, f);
                r.c = 0;//ON
                intersect.push_back(r);
                m.Vertices.push_back(r);
                new_face.Indices.push_back(index);
                
                PushIndex(m.Faces[i], index, j);// вставка значения index на место j
                index++;
            }
        }

        DeleteIndexes(m, m.Faces[i], -1);
        if(out == 4)
            DeleteFace(m, m.Faces[i]);
    }
    DeleteDuplicateVertecies(m, new_face.Indices);
    
    intersect = tries(intersect, f);//вектор вершин в нужном порядке для новой грани
    for(int i = 0; i<intersect.size(); ++i){
        new_face.Indices.push_back(getVertexIndex(intersect[i], m));
    }
    
    m.Faces.push_back(new_face);

    return m;
}

void Triangulation(Mesh &m){
    for(int i = 0; i < m.Faces.size(); ++i){
        if(m.Faces[i].Indices.size() == 4){ // тк модель триангулированная изначально, то кроме сечения максимум 4 индекса у грани
            Face new_face;
            for(int j = 0; j < 3; ++j){
                new_face.Indices.push_back(m.Faces[i].Indices[j]);
                m.Faces[i].Indices.erase(m.Faces[i].Indices.begin() + 1);
            }
            m.Faces.push_back(new_face);
        }
    }
}