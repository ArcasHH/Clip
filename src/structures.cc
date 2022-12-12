#include "structures.h"
#include "obj_parser.h"

#include <unordered_set>
#include <cmath>
#include <algorithm>


double cos(Vector const &v1, Vector const &v2){
    return (v1.dot(v2)) / (v1.length() * v2.length());
}

double PointInFlat (Vertex const &p, Flat const &f){//подставление координат точки в ур.плоскости
    return f.n.x * p.x + f.n.y * p.y + f.n.z * p.z + f.D();
}

void  PointClassify(Mesh &m , Flat const &f){//классификация точек тела относительно плоскости
    double e = 1e-5;
    for(int i =0; i<m.Vertices.size(); ++i){
        if( PointInFlat(m.Vertices[i], f) < -e)//IN
            m.Vertices[i].c = 1;
        else if(PointInFlat(m.Vertices[i], f) < e && PointInFlat(m.Vertices[i], f)  > -e)//ON
            m.Vertices[i].c = 0;
        else if( PointInFlat(m.Vertices[i], f) > e)//OUT
            m.Vertices[i].c = -1;
    }
}

int InPoints(Mesh const &m , Flat const &f){//Количество точек in
    int in = 0;
    for(int i =0; i<m.Vertices.size(); ++i)
        if( m.Vertices[i].c == 1)//IN
            in ++;
    return in;
}
int OnPoints(Mesh const &m , Flat const &f){//Количество точек on 
    int zero = 0;
    for(int i =0; i<m.Vertices.size(); ++i)
        if(  m.Vertices[i].c == 0)//ON
            zero ++;
    return zero;
}
int OutPoints(Mesh const &m , Flat const &f){//Количество точек out
    int out = 0;
    for(int i =0; i<m.Vertices.size(); ++i)
        if(  m.Vertices[i].c == -1)//OUT
            out ++;
    return out;
}

Vertex Segment_Flat_Intersection(Segment const &s, Flat const &f){//точка пересечения плоскости и отрезка
    if(PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0) {
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


void DeleteUncorrectFaces(Mesh &m){
    std::vector<int> index;           //удаление граней, где меньше трех вершин.
    for(int i = 0; i < m.Faces.size() ; ++i)
        if(m.Faces[i].Indices.size() < 3)
            index.push_back(i);

    for(int i = 0; i < index.size(); ++i)
        m.Faces.erase(m.Faces.begin() + index[i] - i);
}
            

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat const &f){//упорядочивание вектора вершин в порядке обхода сечения(соотв.f)

    double e = 1e-9;

    std::vector<Vertex> tries;
    Vector norm = f.n.normalize();
    Vertex geom = {0, 0, 0};

    // std::cout<<"intersect   ";
    // for(int i =0; i < intersect.size(); ++i){
    //     std::cout<<intersect[i].x<<' '<<intersect[i].y<<' '<< intersect[i].z<<' ';
    // }
    // std::cout<<std::endl;

    std::vector<double> angle;

    for(int i=0; i<intersect.size(); ++i){
        geom.x += intersect[i].x / intersect.size();
        geom.y += intersect[i].y / intersect.size();
        geom.z += intersect[i].z / intersect.size();
    }
        //std::cout<< "geom:      "<<geom.x<<' '<<geom.y<<' '<< geom.z<<std::endl;
    double c;
    Vector zero = Vector{geom, intersect[0]};
    for( int i = 1; i < intersect.size(); ++i ){
        c = cos(zero, Vector{geom, intersect[i]});
        //std::cout<<"cos   "<< i<< " cos:    "<<c<<std::endl;;
        Vector Try = zero.cross(Vector{geom, intersect[i]});
        if(Try.length_sq() > e){
            Try = Try.normalize();
            if (Try == norm)
                angle.push_back(c);
            else if (Try == -norm)
                angle.push_back( - 2 - c) ; 
        }
        else
            angle.push_back( -1 );
    }

    // std::cout<<"angle:    ";
    // for (int i =0; i < angle.size(); ++i){
    //     std::cout<<' '<< angle[i]<<' ';
    // }
    // std::cout<<std::endl;

    tries.push_back(intersect[0]);
    sort(angle.begin(), angle.end());
    std::vector<double> angle_sorted;
    for(int i = angle.size()-1; i >=0; i--)
        angle_sorted.push_back(angle[i]);

    angle = angle_sorted;

    //     std::cout<<"angle:    ";
    // for (int i =0; i < angle.size(); ++i){
    //     std::cout<<' '<< angle[i]<<' ';
    // }
    // std::cout<<std::endl;

    for( int i = 0; i < angle.size(); ++i )
        for(int j =0; j < intersect.size(); ++j){
            c = cos(zero, Vector{geom, intersect[j]});
            Vector Try = zero.cross(Vector{geom, intersect[j]});
            if(Try.length_sq() > e){
                Try = Try.normalize();
                if (Try == norm && c == angle[i])
                    tries.push_back(intersect[j]);   
                else if (Try == -norm && (- 2 - c) == angle[i])
                    tries.push_back(intersect[j]);
            }
            else if(Try.length_sq() <= e && std::abs(c-angle[i]) < e)
                tries.push_back(intersect[j]);
        }
    
    return tries;
}


bool SpecialCases(Mesh &m, Flat const &f){
    int zero = OnPoints(m, f);
    int sum = InPoints(m, f) - OutPoints(m, f);

    std::cout<<"in:   "<< InPoints(m,f)<<std::endl;
    std::cout<<"on:   "<< OnPoints(m,f)<<std::endl;
    std::cout<<"out:   "<< OutPoints(m,f)<<std::endl;
    std::cout<<"sum:   "<< sum<<std::endl;

    if(sum == - (m.Vertices.size() - zero) ){//случаи, где результат - пустое множество
        std::cerr << "Empty intersect " << std::endl;
        int size = m.Faces.size();
        // for(int i =0; i < size; ++i)
        //     m.Faces.erase(m.Faces.begin());
        // return true;
    }
    if((sum == (m.Vertices.size() - zero)) ){// случаи, когда результат - исходный объект
        std::cout << "the object will not change " << std::endl;
        return false;
    }
    else
        return true;
}


Mesh ResultOfIntersect( Mesh const &m_in, Flat const &f){
    Mesh m{m_in};
    int index = m.Vertices.size();//индексы новых вершин
    Face new_face;
    std::vector<int> face_del;
    std::vector<Vertex> intersect; //points of intersect

    for(int i = 0; i < m.Faces.size(); ++i){
        int out = 0;
        for(int j =0; j < m.Faces[i].Indices.size(); ++j){
            if(m.Vertices[m.Faces[i].Indices[j]].c == -1)
                out++;

            if( m.Vertices[m.Faces[i].Indices[j]].c  && m.Vertices[m.Faces[i].Indices[j]].c == -m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Segment s;
                s.A = m.Vertices[m.Faces[i].Indices[j]];
                s.B = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r = Segment_Flat_Intersection(s, f);
                r.c = 0;//ON
                m.Vertices.push_back(r);
                m.Faces[i].Indices.insert( std::next(m.Faces[i].Indices.begin(), j + 1) , index);//вставка индекса
                index++;
            }
        }
        DeleteIndexes(m, m.Faces[i], -1);

        if(out == m.Faces[i].Indices.size() + 1)
            face_del.push_back(i);
    }

 



    int del = 0;
    for(int i = 0; i < m.Faces.size(); ++i)//удаление граней из меша, которые содержат 0 вершин(поидее таких не должно возникать)
        if(m.Faces[i].Indices.size() == 0){
            DeleteFace(m, m.Faces[i - del]);
            del++;
        }
                


    std::vector<double> vert;//vector всех координат до удаления из меша чтобы восстаановить индексы для граней

    for(int i = 0; i < m.Vertices.size(); ++i){
        vert.push_back(m.Vertices[i].x);
        vert.push_back(m.Vertices[i].y);
        vert.push_back(m.Vertices[i].z);
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
            for(int i = 0; i < vert.size(); i += 3)
                if(m.Faces[j].Indices[n] == i/3){
                    Vertex v = {{vert[i]}, {vert[i + 1]}, {vert[i + 2]}};
                    m.Faces[j].Indices[n] = getVertexIndex(v, m);
                    break;
                }

    for(int i =0; i < m.Vertices.size(); ++i)//индексы сечения
        if(m.Vertices[i].c == 0)
            new_face.Indices.push_back(i);

    for (int i = 0; i < new_face.Indices.size(); ++i)
        intersect.push_back(m.Vertices[new_face.Indices[i]]);



    intersect = tries( intersect, f);//вектор вершин в нужном порядке для новой грани


    Face face_intersect;//новая грань пересечения
    for(int i =0; i < intersect.size(); ++i)
        face_intersect.Indices.push_back(getVertexIndex(intersect[i], m));

    m.Faces.push_back(face_intersect);
//    for(int i =0; i < m.Vertices.size(); ++i){
//         std::cout<<"i  "<<i<<"     "<<m.Vertices[i].x<<' '<< m.Vertices[i].y<<' '<< m.Vertices[i].z<<"   c:   "<<m.Vertices[i].c<<std::endl;
//     }
//     for(int i =0; i < m.Faces.size(); ++i){
//         for(int j = 0; j < m.Faces[i].Indices.size(); ++j){
//             std::cout<<m.Faces[i].Indices[j]<<' ';
//         }
//         std::cout<<std::endl;
//     }
    return m;

}

void Triangulation(Mesh &m) {



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


bool Check(Mesh const &m){
    std::cout<<"v: "<<m.Vertices.size() <<" f:  "<<m.Faces.size()<<"  res:  "<< m.Vertices.size() - m.Faces.size()/2<<std::endl;
    if( m.Vertices.size() - m.Faces.size()/2 == 2)
        return true;
    else
        return false;
}

void Intersect(Mesh &m, Flat const &f){
    PointClassify(m, f);
    if(SpecialCases(m, f)){
        Mesh res = ResultOfIntersect(m, f);
        Triangulation(res);
        m = res;
    }
    //Check(m);
}