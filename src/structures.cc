#include "structures.h"
#include "obj_parser.h"

#include <learnopengl/model.h>

#include <unordered_set>
#include <cmath>
#include <algorithm>

Mesh Convert(mdl::Mesh const &mesh){//преобразование меша из mdl 
    Mesh m;
    for(int i = 0; i  <mesh.vertices.size(); ++i){
        Vertex v;
        v.x = mesh.vertices[i].Position.x;
        v.y = mesh.vertices[i].Position.y;
        v.z = mesh.vertices[i].Position.z;
        m.Vertices.push_back(v);
    }
    for(int i = 0; i < mesh.indices.size() ; i +=3 ){
        Face new_face;
        for(int n = 0; n < 3; ++n)
            new_face.Indices.push_back(mesh.indices[i + n]) ;
        m.Faces.push_back(new_face);
    }
    return m;
}

double cos(Vector const &v1, Vector const &v2){
    return (v1.dot(v2)) / (v1.length() * v2.length());
}

Flat FlatByPoints(Vertex const&v0, Vertex const&v1, Vertex const&v2){ // построение плоскости по трем точкам
    Flat f;
    f.p = v0;
    f.n.x = -(v1.y - v0.y) * (v2.z - v0.z) + (v1.z - v0.z) * (v2.y - v0.y);
    f.n.y = -(v2.x - v0.x) * (v1.z - v0.z) + (v2.z - v0.z) * (v1.x - v0.x);
    f.n.z = -(v1.x - v0.x) * (v2.y - v0.y) + (v1.y - v0.y) * (v2.x - v0.x);
    f.n = f.n.normalize();
    return f;
}

Vector Normal( Face const &f, Mesh const &m, double precise){ //получение нормали грани(еще можно обобщить)
    if(f.Indices.size() >= 3){
        Vector n = Vector{{m.Vertices[f.Indices[0]]},{m.Vertices[f.Indices[1]]}}.cross(Vector{{m.Vertices[f.Indices[1]]},{m.Vertices[f.Indices[2]]}});
        if( n.length() == 0)
            n = {0, 0, 1};//случаи, когда грань - отрезок
        return n.normalize();
    }
    else
        std::cerr<<"invalid input"<<std::endl;
}

double PointInFlat (Vertex const &p, Flat const &f){//подставление координат точки в ур.плоскости
    return f.n.x * p.x + f.n.y * p.y + f.n.z * p.z + f.D();
}

void  PointClassify(Mesh &m , Flat const &f, double precise){//классификация точек тела относительно плоскости с заданной точностью
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < -precise)//IN
            m.Vertices[i].c = 1;
        else if(PointInFlat(m.Vertices[i], f) < precise && PointInFlat(m.Vertices[i], f)  > -precise)//ON
            m.Vertices[i].c = 0;
        else if( PointInFlat(m.Vertices[i], f) > precise)//OUT
            m.Vertices[i].c = -1;
    }
}

bool SpecialCases(Mesh &m){
    int in_points = 0;
    int on_points = 0;
    int out_points = 0;

    for(int i =0; i<m.Vertices.size(); ++i)
        if( m.Vertices[i].c == 1)//IN
            in_points ++;
        else if( m.Vertices[i].c == 0)//ON
            on_points ++;
        else if( m.Vertices[i].c == -1)//OUT
            out_points ++;
    int sum = in_points - out_points;
    if(sum == - (m.Vertices.size() - on_points) ){//случаи, где результат - пустое множество
        std::cerr << "Empty intersect " << std::endl;
        DeleteMesh( m );
        return false;
    }
    if((sum == (m.Vertices.size() - on_points)) ){// случаи, когда результат - исходный объект
        std::cout << "the object will not change " << std::endl;
        return false;
    }
    else
        return true;
}

Vertex Segment_Flat_Intersection(Segment const &s, Flat const &f){//точка пересечения плоскости и отрезка
    if(PointInFlat(s.A, f) * PointInFlat(s.B, f) < 0) {//если точки по разные стороны относительно плоскости
        Vertex p = s.A - s.B; 
        if((PointInFlat(p, f) - f.D()) != 0){
            double t = - ( PointInFlat(s.B, f)) / (PointInFlat(p, f) - f.D());
            p.x = (s.A.x * t) + (s.B.x * (1.f - t));
            p.y = (s.A.y * t) + (s.B.y * (1.f - t));
            p.z = (s.A.z * t) + (s.B.z * (1.f - t));
            return p;
        }
    }
    throw std::runtime_error("Here should be intersection!!!");  //PointInFlat(p, f) - f.D() == 0 <=> A == B, но тогда невозможно, что PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0
}                                                                // эта функция вызывается для пар точек IN и OUT => PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0 всегда выполнено

int getVertexIndex(Vertex const &v, Mesh const &m){
    for(int i =0; i < m.Vertices.size(); ++i)
        if(m.Vertices[i] == v)
            return i;
    return -1;
}

void DeleteMesh(Mesh &m ){
    m.Faces.clear();
    m.Vertices.clear();
    // int size = m.Faces.size();
    // for(int i =0; i < size; ++i)
    //     m.Faces.erase(m.Faces.begin());
    // size = m.Vertices.size();
    // for(int i =0; i < size; ++i)
    //     m.Vertices.erase(m.Vertices.begin());
}

void DeleteVertex(Mesh &m, Vertex &v){//удаление из массива 1 вершины, совпадающей с v
    std::vector<int> index_to_del;
    for(int i = 0; i < m.Vertices.size(); ++i)
        if( m.Vertices[i] == v){
            m.Vertices.erase(m.Vertices.begin() + i);
            return;
        }
}

void DeleteFace(Mesh &m, Face &f){  //удаление 1 грани из списка
    std::vector<int> index_to_del;
    for(int i = 0; i<m.Faces.size(); i++)
        if (m.Faces[i] == f)
            index_to_del.push_back(i);
    for(int j = 0; j < index_to_del.size(); ++j)
        m.Faces.erase(m.Faces.begin() + index_to_del[j] - j);
}

void DeleteIndexes(Mesh &m, Face &f, int code){//удаление из индексов грани всех, соответсвующих вершинам OUT
    std::vector<int> index_to_del;
    for(int i = 0; i < f.Indices.size(); ++i)
        if( m.Vertices[f.Indices[i]].c == code)
            index_to_del.push_back(i);
    for(int j = 0; j < index_to_del.size(); ++j)
        f.Indices.erase(f.Indices.begin() + index_to_del[j] - j); // проверить еще
}

void DeleteDuplicates(Mesh &m){//удаление дупликатов вершин из меша
    std::vector<int> duplicates;
    std::vector<int> copy;
    for(int i = 0; i < m.Vertices.size() - 1; ++i)
        for(int j = i + 1; j < m.Vertices.size(); ++j)
            if(m.Vertices[i] == m.Vertices[j])
                duplicates.push_back(j);   
    for(int i = 0; i < duplicates.size(); ++i) {
        bool a = true;
        for(int j =0; j < copy.size(); ++j)
            if(copy[j] == duplicates[i])
                a = false;
        if(a)
            copy.push_back(duplicates[i]);
    }
    sort(copy.begin(), copy.end());
    duplicates = copy;
    for(int i = 0; i < duplicates.size(); ++i)
        DeleteVertex(m, m.Vertices[duplicates[i] - i]);
}

void DeleteUncorrectFaces(Mesh &m){//удаление граней, где меньше трех вершин.
    std::vector<int> index;           
    for(int i = 0; i < m.Faces.size() ; ++i)
        if(m.Faces[i].Indices.size() < 3)
            index.push_back(i);
    for(int i = 0; i < index.size(); ++i)
        m.Faces.erase(m.Faces.begin() + index[i] - i);
}
            
std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat const &f, double precise){//упорядочивание вектора вершин в порядке обхода сечения(соотв.f)

    precise = precise * precise;
    std::vector<Vertex> tries;
    Vector norm = f.n.normalize();

    Vertex geom = {0, 0, 0};
    double size = intersect.size();
    for(int i=0; i< size; ++i){
        geom.x += intersect[i].x / size;
        geom.y += intersect[i].y / size;
        geom.z += intersect[i].z / size;
    }

    std::vector<double> angle;
    double c;
    Vector zero = Vector{geom, intersect[0]};
    for( int i = 1; i < intersect.size(); ++i ){//заполнение angle значениями, определяющих положение относительно [0]
        c = cos(zero, Vector{geom, intersect[i]});
        Vector Try = zero.cross(Vector{geom, intersect[i]});
        if(Try.length() > precise){
            Try = Try.normalize();
            if (Try == norm)
                angle.push_back(c);
            else if (Try == -norm)
                angle.push_back( - 2 - c) ; 
        }
        else
            angle.push_back( -1 );
    }

    tries.push_back(intersect[0]);

    sort(angle.begin(), angle.end());//сортировка по углу относительно geom[0]
    std::vector<double> angle_sorted;
    for(int i = angle.size()-1; i >=0; i--)//в обратном порядке
        angle_sorted.push_back(angle[i]);
    angle = angle_sorted;

    for( int i = 0; i < angle.size(); ++i ) //заполнение tries в соответствие с вектором angle (в порядке обхода грани)
        for(int j =0; j < intersect.size(); ++j){
            c = cos(zero, Vector{geom, intersect[j]});
            Vector Try = zero.cross(Vector{geom, intersect[j]});
            if(Try.length() > precise){
                Try = Try.normalize();
                if (Try == norm && std::abs(c-angle[i]) < precise)
                    tries.push_back(intersect[j]);   
                else if (Try == -norm && std::abs( -2 - c - angle[i]) < precise)
                    tries.push_back(intersect[j]);
            }
            else if(Try.length() <= precise && std::abs(c-angle[i]) < precise)
                tries.push_back(intersect[j]);
        }
    return tries;
}

Mesh ResultOfIntersect( Mesh const &m_in, Flat const &f, double precise){
    Mesh m{m_in};
    int index = m.Vertices.size();//индексы новых вершин
    Face new_face;
    std::vector<int> face_del;
    std::vector<Vertex> intersect; //points of intersect

    for(int i = 0; i < m.Faces.size(); ++i){
        for(int j =0; j < m.Faces[i].Indices.size(); ++j){//для каждой грани обходим все индексы

            if( m.Vertices[m.Faces[i].Indices[j]].c  && m.Vertices[m.Faces[i].Indices[j]].c == -m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Segment s;                      //если вершина j не ON, и следующая за j(IN) - OUT или наоборот, то строим пересечение ребра и плоскости
                s.A = m.Vertices[m.Faces[i].Indices[j]];
                s.B = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r = Segment_Flat_Intersection(s, f);// то строим пересечение ребра и плоскости
                r.c = 0;//ON
                m.Vertices.push_back(r);//добавляем новую вершину в список вершин
                m.Faces[i].Indices.insert( std::next(m.Faces[i].Indices.begin(), j + 1) , index);//вставка индекса новой вершины в список индексов грани
                index++;
            }
        }
        DeleteIndexes(m, m.Faces[i], -1);//удаление из граней индексов, соответствующим вершинам OUT
    }

    DeleteUncorrectFaces(m);

   //vector всех вершин до удаления из меша чтобы восстаановить индексы для граней
    std::vector<Vertex> vert;
    for(int i = 0; i < m.Vertices.size(); ++i){
        vert.push_back(m.Vertices[i]);
    }

    DeleteDuplicates(m);//удаление совпадающих вершин

    std::vector<int> index_ver_del;//удаление вершин OUT из меша
    for(int i = 0; i < m.Vertices.size(); ++i){
        if(m.Vertices[i].c == -1){
            index_ver_del.push_back(i);
        }
    }
    for(int j = 0; j < index_ver_del.size(); ++j){
        DeleteVertex(m, m.Vertices[index_ver_del[j] -j]);
    }

    for(int j =0; j < m.Faces.size(); ++j) //переопределение индексов вершин для граней после удаления вершин
        for(int n =0; n < m.Faces[j].Indices.size(); ++n)
            for(int i = 0; i < vert.size(); i ++)
                if(m.Faces[j].Indices[n] == i){
                    m.Faces[j].Indices[n] = getVertexIndex(vert[i], m);
                    break;
                }

    for(int i =0; i < m.Vertices.size(); ++i)//индексы сечения
        if(m.Vertices[i].c == 0)
            new_face.Indices.push_back(i);

    for (int i = 0; i < new_face.Indices.size(); ++i)
        intersect.push_back(m.Vertices[new_face.Indices[i]]);

    intersect = tries( intersect, f, precise);//вектор вершин intersect в нужном порядке для новой грани

    Face face_intersect;//новая грань пересечения/ заполнение индексами вершин intersect
    for(int i =0; i < intersect.size(); ++i)
        face_intersect.Indices.push_back(getVertexIndex(intersect[i], m));

    m.Faces.push_back(face_intersect);

    return m;

}


void Triangulation(Mesh &m) { //переделай на триангуляцию каждой грани

    DeleteUncorrectFaces(m);

    int faces = m.Faces.size();
    for(int i = 0; i < faces ; ++i){
        if(m.Faces[i].Indices.size() == 3) continue;
        std::vector<int> index;
        for(int j = 0; j < m.Faces[i].Indices.size() - 1; j += 2){
            Face new_face;
            for( int n = 0; n < 3; ++n)
                new_face.Indices.push_back(m.Faces[i].Indices[(j + n)%m.Faces[i].Indices.size()]);
            index.push_back(j+1);
            m.Faces.push_back(new_face);
        }
        int size = m.Faces[i].Indices.size();
        if(size > 3)
            for(int j = size/2 ; j > 0; j-- )
                m.Faces[i].Indices.erase(m.Faces[i].Indices.begin() + j + j - 1);
    }

    for(int i =0; i < m.Faces.size(); ++i){ ///удаление внутри циклв\а
        for(int j =0; j < m.Faces[i].Indices.size() - 1; ++j)
            for(int k = j+1; k < m.Faces[i].Indices.size(); ++k)
                if(m.Faces[i].Indices[j] == m.Faces[i].Indices[k])
                    m.Faces[i].Indices.erase(m.Faces[i].Indices.begin() + k);      
    }
    DeleteUncorrectFaces(m);
    for(int i = 0; i < m.Faces.size(); ++i)//если остались грани размера больше 3, то вызывается рекурсивно
        if(m.Faces[i].Indices.size() > 3){
            Triangulation(m);
            break;
        }
}

bool Check(Mesh const &m){//проверка н аправильность построения формулой Эйлера
        return (m.Vertices.size() - m.Faces.size()/2 == 2);
}

void Intersect(Mesh &m, Flat const &f, double precise){ //отсечение объекта плоскостью
    PointClassify(m, f, precise);
    if(SpecialCases(m)){
        Mesh res = ResultOfIntersect(m, f, precise);
        //Correct(res, precise);
        Triangulation(res);
        m = res;
    }
}

void Correct(Mesh &m, double precise){
    // DeleteUncorrectFaces(m);
    // std::cout<<std::endl;
    // std::vector<Vector> normals;
    // Mesh m2{m};
    
    // for(int i =0; i < m.Faces.size(); ++i){
    //     Vector n = Vector{{m.Vertices[m.Faces[i].Indices[0]]},{m.Vertices[m.Faces[i].Indices[1]]}}.cross(Vector{{m.Vertices[m.Faces[i].Indices[1]]},{m.Vertices[m.Faces[i].Indices[2]]}});
    //     bool uniq = true;
    //     n = n.normalize();
    //     if(normals.size() != 0)
    //         for(int j =0; j < normals.size(); ++j)
    //             if(normals[j] == n)
    //                 uniq = false; 
    //     if(uniq)
    //         normals.push_back(n);
    // }

    // std::vector<Face> faces;
    // for(int i =0; i < normals.size(); ++i){
    //     Face new_face;
    //     std::vector<Vertex> vert;

    //     for(int j =0; j < m.Faces.size(); ++j){
    //         //std::cout<<"faces   "<<m.Vertices[m.Faces[j].Indices[0]];
    //         Vector n = Vector{{m.Vertices[m.Faces[j].Indices[0]]},{m.Vertices[m.Faces[j].Indices[1]]}}.cross(Vector{{m.Vertices[m.Faces[j].Indices[1]]},{m.Vertices[m.Faces[j].Indices[2]]}});
    //         if(n.length() >= precise)
    //             n = n.normalize();
    //         //std::cout<<"  n  "<< n <<"   normal   "<<normals[i]<<std::endl;
    //         if(n == normals[i]){
    //             for(int k =0; k < m.Faces[j].Indices.size(); ++k){
    //                 bool uniq = true;
    //                 for(int n =0; n < vert.size(); ++n)
    //                     if(vert[n] == m.Vertices[m.Faces[j].Indices[k]])
    //                         uniq = false;
    //                 if(uniq)
    //                     vert.push_back(m.Vertices[m.Faces[j].Indices[k]]);
    //             }
    //         }
    //     }

    //     Flat f;
    //     f.n = normals[i];
    //     f.p = {0, 0, 0};

    //     vert = tries( vert, f, precise);
    
    //     int size = vert.size();
    //     for(int k =0; k < size ; ++k){
    //         Vector n = Vector{{vert[k % vert.size()]},{vert[(k+1) % vert.size()]}}.cross(Vector{vert[(k+1) % vert.size()],vert[(k+2) % vert.size()]});
    //         if(n.length_sq() < precise){
    //             DeleteVertex(m2, vert[(k+1) % vert.size()]);
    //             vert.erase(vert.begin() + (k + 1) % vert.size());
    //         }
    //     }

    //     for(int j =0; j < vert.size(); ++j)
    //         new_face.Indices.push_back(getVertexIndex(vert[j], m));
    //     faces.push_back(new_face);
    // }

    // for(int i =0; i < faces.size(); ++i)
    //     for(int j = 0; j < faces[i].Indices.size(); ++j)
    //         faces[i].Indices[j] = getVertexIndex( m.Vertices[faces[i].Indices[j]] , m2);

    // m.Faces = faces;
    // m.Vertices = m2.Vertices;
}



void ClassifyObjects(Mesh &m1, Mesh &m2, double precise){//Классификация точек m1 относительно m2
    // for(int i =0; i<m1.Vertices.size(); ++i){
    //     int in_points = 0;
    //     int on_points = 0;
    //     int out_points = 0;
    //     for(int j =0; j < m2.Faces.size(); ++j){
    //         Flat f = FlatByPoints(m2.Vertices[m2.Faces[j].Indices[0]], m2.Vertices[m2.Faces[j].Indices[2]], m2.Vertices[m2.Faces[j].Indices[1]]);
    //         if( PointInFlat(m1.Vertices[i], f) < -precise)//IN
    //             in_points += 1;
    //         else if(PointInFlat(m1.Vertices[i], f) < precise && PointInFlat(m1.Vertices[i], f)  > -precise){//ON
    //             on_points += 1;
    //         }
    //         else if( PointInFlat(m1.Vertices[i], f) > precise){//OUT
    //             out_points += 1;
    //         }
    //     }
    //     std::cout<<"points:   "<<in_points<<' '<<on_points<<' '<<out_points<<std::endl;
    //     if(in_points == m2.Faces.size() && out_points == 0){
    //         m1.Vertices[i].c = 1;
    //     }
    //     else if(on_points + in_points == m2.Faces.size() && out_points == 0){
    //         m1.Vertices[i].c = 0;
    //     }
    //     else if(out_points != 0){
    //         m1.Vertices[i].c = -1;
    //     }
    // }
}


Mesh ResultOfDifference(Mesh &m1, Mesh &m2, double precise){

}

void bool_union(Mesh &m1, Mesh &m2, double precise){
    // ClassifyObjects(m1, m2, precise);

    // std::cout<< "m1:  "<<std::endl;
    // for(int i=0; i < m1.Vertices.size(); ++i){
    //     std::cout<< m1.Vertices[i].x << ' '<<m1.Vertices[i].y << ' '<<m1.Vertices[i].z << "   c:    "<<m1.Vertices[i].c << ' '<<std::endl;
    // }
    // std::cout<< "m2:  "<<std::endl;
    //     for(int i=0; i < m2.Vertices.size(); ++i){
    //     std::cout<< m2.Vertices[i].x << ' '<<m2.Vertices[i].y << ' '<<m2.Vertices[i].z << "   c:    "<<m2.Vertices[i].c << ' '<<std::endl;
    // }

    // if(SpecialCases(m1)){
    //     Mesh res = ResultOfDifference(m1, m2, precise);
    //     Triangulation(res);
    //     m1 = res;
    // }
}

void bool_difference(Mesh &m1, Mesh &m2, double precise){

    // ClassifyObjects(m1, m2, precise); //классификация m1 относительно m2
    // ClassifyObjects(m2, m1, precise);//классификация m2 относительно m1

    // std::cout<< "m1:  "<<std::endl;
    // for(int i=0; i < m1.Vertices.size(); ++i){
    //     std::cout<< m1.Vertices[i].x << ' '<<m1.Vertices[i].y << ' '<<m1.Vertices[i].z << "   c:    "<<m1.Vertices[i].c << ' '<<std::endl;
    // }
    // std::cout<< "m2:  "<<std::endl;
    //     for(int i=0; i < m2.Vertices.size(); ++i){
    //     std::cout<< m2.Vertices[i].x << ' '<<m2.Vertices[i].y << ' '<<m2.Vertices[i].z << "   c:    "<<m2.Vertices[i].c << ' '<<std::endl;
    // }

    // if(SpecialCases(m1)){
    //     Mesh res = ResultOfDifference(m1, m2, precise);
    //     Triangulation(res);
    //     m1 = res;
    // }
}

