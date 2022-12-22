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
    if(f.n.length() < precise)
        f.n = {0, 0, 1};
    f.n = f.n.normalize();
    
    return f;
}

Vector Normal( Face const &f, Mesh const &m, double precise){ //получение нормали грани(еще можно обобщить)
    if(f.Indices.size() >= 3){
        Vector n = Vector{{m.Vertices[f.Indices[0]]},{m.Vertices[f.Indices[1]]}}.cross(Vector{{m.Vertices[f.Indices[1]]},{m.Vertices[f.Indices[2]]}});
        if( n.length() < precise)
            n = {0, 0, 1};//случаи, когда грань - отрезок
        return n.normalize();
    }
    else
        return Vector{0, 0, 0};
        // std::cerr<<"invalid input"<<std::endl;
}

double PointInFlat (Vertex const &p, Flat const &f){//подставление координат точки в ур.плоскости
    return f.n.x * p.x + f.n.y * p.y + f.n.z * p.z + f.D();
}

void  PointClassify(Mesh &m , Flat const &f, double precise){//классификация точек тела относительно плоскости с заданной точностью
    precise = precise * 10;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < -precise)//IN
            m.Vertices[i].c = 1;
        else if(std::abs(PointInFlat(m.Vertices[i], f)) <= precise)//ON
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
            p.x = (s.A.x * t) + (s.B.x * (1 - t));
            p.y = (s.A.y * t) + (s.B.y * (1. - t));
            p.z = (s.A.z * t) + (s.B.z * (1. - t));
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
    for(int i = 0; i < m.Faces.size(); i++)
        if (m.Faces[i] == f){
            index_to_del.push_back(i);
            m.Faces[i].Indices.clear();
        }
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

    std::vector<double> angle;//будем смотреть на угол между вектором, соед. геом центр g и intersect[0], и вектором,  соед. геом центр g и intersect[i], i != 0
    double c;
    Vector zero = Vector{geom, intersect[0]};
    for( int i = 1; i < intersect.size(); ++i ){//заполнение angle
        c = cos(zero, Vector{geom, intersect[i]});
        Vector Try = zero.cross(Vector{geom, intersect[i]});
        if(Try.length() > precise){
            Try = Try.normalize();
            if (Try == norm)//sin > 0
                angle.push_back(c);
            else if (Try == -norm)//sin < 0
                angle.push_back( - 2 - c) ; 
        }
        else // sin = 0
            angle.push_back( -1 );
    }

    tries.push_back(intersect[0]);

    sort(angle.begin(), angle.end(), std::greater<double>());//сортировка по углу

    for( int i = 0; i < angle.size(); ++i ) //заполнение tries в соответствие с вектором angle (в порядке обхода грани)
        for(int j =0; j < intersect.size(); ++j){
            c = cos(zero, Vector{geom, intersect[j]});
            Vector Try = zero.cross(Vector{geom, intersect[j]});
            if(Try.length() > precise){
                Try = Try.normalize();
                if (Try == norm && std::abs(c-angle[i]) < precise)//sin > 0
                    tries.push_back(intersect[j]);   
                else if (Try == -norm && std::abs( -2 - c - angle[i]) < precise)//sin < 0
                    tries.push_back(intersect[j]);
            }
            else if(std::abs(c-angle[i]) < precise)//sin == 0
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
        for(int j =0; j < m.Faces[i].Indices.size(); ++j)//для каждой грани обходим все индексы
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
        DeleteIndexes(m, m.Faces[i], -1);//удаление из граней индексов, соответствующим вершинам OUT
    }

    DeleteUncorrectFaces(m);

   //vector всех вершин до удаления из меша чтобы восстаановить индексы для граней
    std::vector<Vertex> vert;
    for(int i = 0; i < m.Vertices.size(); ++i)
        vert.push_back(m.Vertices[i]);

    DeleteDuplicates(m);//удаление совпадающих вершин

    std::vector<int> index_ver_del;//удаление вершин OUT из меша
    for(int i = 0; i < m.Vertices.size(); ++i)
        if(m.Vertices[i].c == -1)
            index_ver_del.push_back(i);
    
    for(int j = 0; j < index_ver_del.size(); ++j)
        DeleteVertex(m, m.Vertices[index_ver_del[j] -j]);

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
    for (int i = 0; i < new_face.Indices.size(); ++i)//заполнение intersect
        intersect.push_back(m.Vertices[new_face.Indices[i]]);

    intersect = tries( intersect, f, precise);//вектор вершин intersect в нужном порядке для новой грани

    Face face_intersect;//новая грань пересечения/ заполнение индексами вершин intersect
    for(int i =0; i < intersect.size(); ++i)
        face_intersect.Indices.push_back(getVertexIndex(intersect[i], m));

    m.Faces.push_back(face_intersect);
    return m;
}


void Triangulation(Mesh &m) { //переделай на триангуляцию каждой грани

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
        if(m.Faces[i].Indices.size() > 3)
            i --;
    }

    for(int i =0; i < m.Faces.size(); ++i) ///удаление дубликатов индексов в грани(если присутствуют по одной копии)
        for(int j =0; j < m.Faces[i].Indices.size() - 1; ++j)
            for(int k = j+1; k < m.Faces[i].Indices.size(); ++k)
                if(m.Faces[i].Indices[j] == m.Faces[i].Indices[k])
                    m.Faces[i].Indices.erase(m.Faces[i].Indices.begin() + k);      

    DeleteUncorrectFaces(m);
}

void Correct(Mesh &m, double precise){

    std::vector<Flat> flats;
    for(int i =0; i < m.Faces.size(); ++i){
        if(Vector{{m.Vertices[m.Faces[i].Indices[0]]},{m.Vertices[m.Faces[i].Indices[1]]}}.cross(Vector{{m.Vertices[m.Faces[i].Indices[1]]},{m.Vertices[m.Faces[i].Indices[2]]}}).length_sq() <= precise){
            for(int n =0; n < 3; ++n){
                m.Faces[i].Indices.erase(m.Faces[i].Indices.begin());
            }
        }
        Flat f;
        f.n = Normal(m.Faces[i], m, precise);
        f.p = m.Vertices[m.Faces[i].Indices[0]];
        m.Faces[i].norm = f.n;
        bool uniq = true;
        for(int j =0; j < flats.size(); ++j){
            if( f == flats[j])
                uniq = false;
        }
        if(uniq && !(f.n == Vector{0, 0, 0})){
            flats.push_back(f);
        }
    }
    DeleteUncorrectFaces(m);

    std::vector< std::vector<Vertex>> vertices;
    for(int i =0; i < flats.size(); ++i){
        std::vector<Vertex> vert;
        for(int j = 0; j < m.Vertices.size(); ++j){
            if(std::abs(PointInFlat(m.Vertices[j], flats[i])) < precise){
                vert.push_back(m.Vertices[j]);
            }
        }
        vert = tries(vert, flats[i], precise);

        int size = vert.size();
        std::vector<Vertex> del;
        std::vector<Vertex> vert_copy;
        for(int j = 0; j < size; ++j){
            Vector n = Vector{{vert[j % vert.size()]},{vert[(j + 1)% vert.size()]}}.cross(Vector{{vert[(j + 1)% vert.size()]},{vert[(j + 2)% vert.size()]}});
            if(n.length() < precise){
                del.push_back(vert[(j + 1)% vert.size()]);
                m.Vertices.erase(m.Vertices.begin() + getVertexIndex(vert[(j + 1)% vert.size()], m));
                //std::cout<<"del     "<< i <<"    "<<vert[(j + 1)% vert.size()]<<std::endl;
            }
        }
        for(int j = 0; j < vert.size(); ++j){
            bool not_deleted = true;
            for(int k = 0; k < del.size(); ++k)
                if(vert[j] == del[k])
                    not_deleted = false;
            if(not_deleted)
                vert_copy.push_back(vert[j]);
        }
        vert = vert_copy;
        vertices.push_back(vert);
    }


    m.Faces.clear();
    for(int i =0; i < vertices.size(); ++i){
        Face new_face;
        for(int j =0; j < vertices[i].size(); ++j)
            new_face.Indices.push_back(getVertexIndex(vertices[i][j], m));
        new_face.norm = flats[i].n;
        m.Faces.push_back(new_face);
    }
}

bool Check(Mesh const &m){//проверка н аправильность построения формулой Эйлера
        return (m.Vertices.size() - m.Faces.size()/2 == 2);
}

void Intersect(Mesh &m, Flat const &f, double precise){ //отсечение объекта плоскостью
    PointClassify(m, f, precise);
    if(SpecialCases(m)){
        Mesh res = ResultOfIntersect(m, f, precise);
        Triangulation(res);
        Correct(res, precise);
        Triangulation(res);
        m = res;
    }
}


void Adj(Mesh & m){
    for(int i =0; i < m.Faces.size(); ++i)
        for(int j =0; j < m.Faces[i].Indices.size(); ++j)
            m.Vertices[m.Faces[i].Indices[j]].faces.push_back(i);
}


void ClassifyObjects(Mesh &m1, Mesh &m2, double precise){//Классификация точек m1 относительно m2(для триангулированных)
    
    for(int i =0; i<m1.Vertices.size(); ++i){
        int in_points = 0;
        int on_points = 0;
        int out_points = 0;
        for(int j =0; j < m2.Faces.size(); ++j){
            Flat f = FlatByPoints(m2.Vertices[m2.Faces[j].Indices[0]], m2.Vertices[m2.Faces[j].Indices[2]], m2.Vertices[m2.Faces[j].Indices[1]]);
            if( PointInFlat(m1.Vertices[i], f) < -precise)//IN
                in_points += 1;
            else if(PointInFlat(m1.Vertices[i], f) < precise && PointInFlat(m1.Vertices[i], f)  > -precise)//ON
                on_points += 1;
            else if( PointInFlat(m1.Vertices[i], f) > precise)//OUT
                out_points += 1;
        }
        if(in_points == m2.Faces.size() && out_points == 0)
            m1.Vertices[i].c = 1;
        else if(on_points + in_points == m2.Faces.size() && out_points == 0)
            m1.Vertices[i].c = 0;
        else if(out_points != 0)
            m1.Vertices[i].c = -1;
    }
}

void bool_intersectiom(Mesh &m1, Mesh &m2, double precise){//пересечение m1 и m2
    ClassifyObjects(m1, m2, precise);
    ClassifyObjects(m2, m1, precise);
    Correct(m1, precise);
    Correct(m2, precise);
    
    Adj( m1 );
    Adj( m2 );

    int in1 =0;
    int in2 = 0;
    int on1 =0;
    int on2 = 0;
    for(int i =0; i < m1.Vertices.size(); ++i){
        if(m1.Vertices[i].c == 1)
           in1 ++;
        if(m1.Vertices[i].c == 0)
           on1 ++;
    }
        
    for(int i =0; i < m2.Vertices.size(); ++i){
        if(m2.Vertices[i].c == 1)
            in2 ++;
        if(m2.Vertices[i].c == 0)
           on2 ++;
    }  
    std::cout<<"  in1:  "<<in1<<"   in2:  "<<in2<<"  on1:  "<<on1<<"   on2:  "<<on2<<std::endl;

    // for(int i =0; i < m1.Vertices.size(); ++i){
    //     std::cout<<"vert 1   "<< i <<"   "<<m1.Vertices[i]<<"     c:    "<<m1.Vertices[i].c<<::endl;
    //     std::cout<<"f1     ";
    //     for(int j =0; j < m1.Vertices[i].faces.size(); ++j){
    //         std::cout<< m1.Vertices[i].faces[j]<<' ';
    //     }
    //     std::cout<<std::endl;
    // }
    // for(int i =0; i < m1.Faces.size(); ++i){
    //     std::cout<<" face 1  "<<i<<"    ";
    //     for(int j =0; j < m1.Faces[i].Indices.size(); ++j){
    //         std::cout<<m1.Faces[i].Indices[j]<<' ';
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;
    // for(int i =0; i < m2.Vertices.size(); ++i){
    //     std::cout<<"vert 2   "<< i <<"   "<<m2.Vertices[i]<<"     c:    "<<m2.Vertices[i].c<<::endl;
    //     std::cout<<"f2    ";
    //     for(int j =0; j < m2.Vertices[i].faces.size(); ++j){
    //         std::cout<< m2.Vertices[i].faces[j]<<' ';
    //     }
    //     std::cout<<std::endl;
    // }
    // for(int i =0; i < m2.Faces.size(); ++i){
    //     std::cout<<" face 2  "<<i<<"    ";
    //     for(int j =0; j < m2.Faces[i].Indices.size(); ++j){
    //         std::cout<<m2.Faces[i].Indices[j]<<' ';
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;

    // std::cout<<"ASD2"<<std::endl;
    for(int i =0; i < m2.Vertices.size(); ++i)
        if(m2.Vertices[i].c != -1){
            for(int j =0; j < m2.Vertices[i].faces.size(); ++j){
                Face f = m2.Faces[m1.Vertices[i].faces[j]];
                Flat flat = FlatByPoints(m2.Vertices[f.Indices[0]], m2.Vertices[f.Indices[2]], m2.Vertices[f.Indices[1]]);
                Intersect(m1, flat, precise);
            }
    }
    Triangulation(m1);
    Triangulation(m2);
}