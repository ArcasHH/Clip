#include "structures.h"
#include "obj_parser.h"

#include <unordered_set>

float PointInFlat (Vertex const &p, Flat const &f){//подставление координат точки в ур.плоскости
    return f.n.x * p.x + f.n.y * p.y + f.n.z * p.z + f.D();
}

void  PointClassify(Mesh &m , Flat const &f){//классификация точек тела относительно плоскости
    float e = 1e-4;
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
            float t = - ( PointInFlat(s.B, f)) / (PointInFlat(p, f) - f.D());
        //p = (s.A * t) + (s.B * (1.f - t));
        p.x = (s.A.x * t) + (s.B.x * (1.f - t));
        p.y = (s.A.y * t) + (s.B.y * (1.f - t));
        p.z = (s.A.z * t) + (s.B.z * (1.f - t));
        return p;
        }
    }
    throw std::runtime_error("Here should be intersection!!!");  //PointInFlat(p, f) - f.D() == 0 <=> A == B, но тогда невозможно, что PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0
}                                                                // эта функция вызывается для пар точек IN и OUT => PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0 всегда выполнено

int getVertexIndex(Vertex const &v, Mesh const &m){ //vector.find
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

void PushIndex(Face &f, int index, int j){//вставка индекса новой вершины в список на место j+1 (после j) в списке грани
    f.Indices.insert( std::next(f.Indices.begin(), j + 1) , index); 
}

// bool isOnLine(Vertex &v, Vertex const &v1, Vertex const &v2){
//     return Vector(v, v1).normalize() == Vector(v2, v).normalize();
// }

// std::vector<int> DuplicateVertecies(Mesh &m, std::vector<int> arr){//удаление по одному из пaр индексов, соответствующих одинаковым вершинам
//     std::vector<int> d;
//     std::vector<int> t;
//     for(int i =0; i < arr.size() - 1; ++i)

//         for(int j = i + 1; j < arr.size(); ++j)
//             if(m.Vertices[arr[i]] == m.Vertices[arr[j]]){
//                 d.push_back(arr[j]);
//                 t.push_back(arr[i]);
//             }     
//     for(int i = 0; i < d.size(); ++i)
//         for(int n = 0; n < m.Faces.size(); ++n)
//             for(int k = 0; k < m.Faces[n].Indices.size(); ++k)
//                 if(m.Faces[n].Indices[k] == d[i])
//                     m.Faces[n].Indices[k] = t[i];
//     return t;
// }
    
void DeleteDuplicates(Mesh &m){//удаление дупликатов вершин из меша
    std::vector<int> duplicates;
    for(int i = 0; i < m.Vertices.size() - 1; ++i)
        for(int j = i + 1; j < m.Vertices.size(); ++j)
            if(m.Vertices[i] == m.Vertices[j])
                duplicates.push_back(j);              
    for(int i = 0; i < duplicates.size(); ++i)
        DeleteVertex(m, m.Vertices[duplicates[i] - i]);
}


void DeleteUncorrectFaces(Mesh &m){
    std::vector<int> index;           //удаление граней, где меньше трех вершин.(почему возникают?)
    for(int i = 0; i < m.Faces.size() ; ++i)
        if(m.Faces[i].Indices.size() < 3)
            index.push_back(i);
    for(int i = 0; i < index.size(); ++i)
        m.Faces.erase(m.Faces.begin() + index[i] - i);
}
            

std::vector<Vertex> tries (std::vector<Vertex> &intersect, Flat const &f){//упорядочивание вектора вершин в порядке обхода сечения(соотв.f)
    std::vector<Vertex> tries;
    Vector norm = f.n.normalize();
    //std::vector<int> del;

    // std::unordered_set<int> unique_id_to_delete;
    // for(int i = 0; i < intersect.size(); ++i)
    //     for(int j = 0; j < intersect.size()-1; ++j)
    //         if( j != i)
    //             for(int k = j+1; k < intersect.size(); ++k)
    //                 if(k != i && k != j)
    //                     if(isOnLine(intersect[i], intersect[j], intersect[k])){
    //                         unique_id_to_delete.insert(i);
    //                         break;
    //                     }
    
    // std::vector<Vertex> new_vertices;
    // for (int i = 0; i != intersect.size(); ++i) {
    //     if (!unique_id_to_delete.count(i))
    //         new_vertices.push_back(intersect[i]);
    // }
    // intersect = new_vertices;

    std::vector<int> arr(intersect.size()-1);//массив - количество появления вершин на второй позиции
    for(int i=0; i<intersect.size() - 1; ++i)
        arr[i] = 0;

    for( int i = 1; i < intersect.size() - 1; ++i )
        for(int j = i + 1; j < intersect.size(); ++j ){
            Vector Try = Vector{intersect[0], intersect[i]}.cross(Vector{intersect[i], intersect[j]}).normalize();
            if (Try == norm)
                arr[i-1]++;
            else if (Try == -norm)
                arr[j-1]++;
        }

    tries.push_back(intersect[0]);
    for(int j = 2; j <= intersect.size(); j++){ //tries - список вершин новой грани в порядке обхода. (сортировка)
        for(int i=0; i < intersect.size() - 1; ++i)
            if(arr[i] == intersect.size() - j)
                tries.push_back(intersect[i+1]);
    }
    return tries;
}





bool SpecialCases(Mesh &m, Flat const &f){//проверка на случаи полного и пустого пересечений
    int zero = OnPoints(m, f);
    int sum = InPoints(m, f) - OutPoints(m, f);

    std::cout<<"in:   "<< InPoints(m,f)<<std::endl;
    std::cout<<"on:   "<< OnPoints(m,f)<<std::endl;
    std::cout<<"out:   "<< OutPoints(m,f)<<std::endl;
    std::cout<<"sum:   "<< sum<<std::endl;

    if(sum == - (m.Vertices.size() - zero) ){//случаи, где результат - пустое множество
        std::cout << "Empty intersect " << std::endl;
        //отрисовать пустой экран
        return false;
    }
    if((sum == (m.Vertices.size() - zero)) ){// случаи, когда результат - исходный объект
        std::cout << "the object will not change " << std::endl;
        //отрисовка исходного меша
        return true;
    }
    else
        return false;
}


Mesh ResultOfIntersect( Mesh const &m_in, Flat const &f){
    Mesh m{m_in};
    int index = m.Vertices.size();//индексы новых вершин
    Face new_face;
    std::vector<int> face_del;
    std::vector<Vertex> intersect; //points of intersect

    // for(int i =0; i < m.Faces.size(); ++i){
    //     std::cout<< "doddd  "<< i <<"   ";
    //     for(int j =0; j < m.Faces[i].Indices.size(); ++j){
    //         std::cout<< m.Faces[i].Indices[j] <<' ';
    //     }
    //     std::cout<<std::endl;
    // }


    // for(int i =0; i<m.Vertices.size(); ++i){
    //     std::cout<< "resssvert  "<< i <<"   "<< m.Vertices[i].x <<' '<< m.Vertices[i].y<<' '<< m.Vertices[i].z<<std::endl;
    // }

    // for(int i =0; i < m.Faces.size(); ++i){
    //     std::cout<< "ressf  "<< i <<"   ";
    //     for(int j =0; j < m.Faces[i].Indices.size(); ++j){
    //         std::cout<< m.Faces[i].Indices[j] <<' ';
    //     }
    //     std::cout<<std::endl;
    // }


    
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
                
                PushIndex(m.Faces[i], index, j);// вставка значения index на место j
                index++;
            }
        }
        DeleteIndexes(m, m.Faces[i], -1);

        if(out == m.Faces[i].Indices.size() + 1)
            face_del.push_back(i);
    }

    
    for(int i = 0; i < face_del.size(); ++i)//удаление граней из меша, все вершины которой OUT
        DeleteFace(m, m.Faces[face_del[i] - i]);
    
    // int del = 0;
    // for(int i = 0; i < m.Faces.size(); ++i){//удаление граней из меша, которые содержат 0 вершин(поидее таких не должно возникать)
    //     if(m.Faces[i].Indices.size() == 0){
    //         DeleteFace(m, m.Faces[i - del]);
    //         del++;
    //     }
    // }

    std::vector<float> vert;//vector всех координат до удаления из меша чтобы восстаановить индексы для граней

    //new_face.Indices = DuplicateVertecies(m, new_face.Indices);
    for(int i = 0; i < m.Vertices.size(); ++i){
        vert.push_back(m.Vertices[i].x);
        vert.push_back(m.Vertices[i].y);
        vert.push_back(m.Vertices[i].z);
    }

    DeleteDuplicates(m);

    std::vector<int> index_ver_del;//удаление вершин OUT из меша
    for(int i = 0; i < m.Vertices.size(); ++i){
        if(m.Vertices[i].c == -1){
            index_ver_del.push_back(i);
        }
    }
    for(int j = 0; j < index_ver_del.size(); ++j){
        DeleteVertex(m, m.Vertices[index_ver_del[j] -j]);
    }


    for(int j =0; j < m.Faces.size(); ++j) 
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



    intersect = tries(intersect, f);//вектор вершин в нужном порядке для новой грани


    Face face_intersect;
    for(int i =0; i < intersect.size(); ++i)
        face_intersect.Indices.push_back(getVertexIndex(intersect[i], m));

    m.Faces.push_back(face_intersect);

    for(int i =0; i < m.Faces.size() - 1; ++i){  /// если несколько дублирующихся граней, то возможны jo
        for(int j =i + 1; j < m.Faces.size(); ++j)
            for(int n =0; n < 3; ++n)
                if(m.Faces[i].Indices[0] == m.Faces[j].Indices[n] && m.Faces[i].Indices[1] == m.Faces[j].Indices[(n+1)% 3 ] && m.Faces[i].Indices[2] == m.Faces[j].Indices[(n + 2)% 3 ])
                    m.Faces.erase(m.Faces.begin() + j);
    }


    return m;

}

void Triangulation(Mesh &m) {


    // for(int i =0; i<m.Vertices.size(); ++i){
    //     std::cout<< "dovert  "<< i <<"   "<< m.Vertices[i].x <<' '<< m.Vertices[i].y<<' '<< m.Vertices[i].z<<std::endl;
    // }

    // for(int i =0; i < m.Faces.size(); ++i){
    //     std::cout<< "dof  "<< i <<"   ";
    //     for(int j =0; j < m.Faces[i].Indices.size(); ++j){
    //         std::cout<< m.Faces[i].Indices[j] <<' ';
    //     }
    //     std::cout<<std::endl;
    // }


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
            // for(int k =0; k < new_face.Indices.size(); k ++){
            //     std::cout<<' '<<new_face.Indices[k]<<' ';
            // }
            // std::cout<<std::endl;
        }
        int size = m.Faces[i].Indices.size();

        if(size > 3)
            for(int j = size/2 ; j > 0; j-- )
                m.Faces[i].Indices.erase(m.Faces[i].Indices.begin() + j + j - 1);
    }

    
    DeleteUncorrectFaces(m);


}

void Correct(Mesh &m, Flat const &f){//перестраивает триангулированный меш в обычный

    Flat plane;
    std::vector<Vertex> vert;
    Face intersect_face;

    // for(int i =0; i < m.Faces.size(); ++i){
    //     std::cout<< "doooplaned 3     "<< i <<"   ";
    //     for(int j =0; j < m.Faces[i].Indices.size(); ++j){
    //         std::cout<< m.Faces[i].Indices[j] <<' ';
    //     }
    //     std::cout<<std::endl;
    // }

    // for(int i =0; i<m.Vertices.size(); ++i){
    //     std::cout<< "verticles in correct  "<< i <<"   "<< m.Vertices[i].x <<' '<< m.Vertices[i].y<<' '<< m.Vertices[i].z<<' '<<m.Vertices[i].c<<std::endl;
    // }

    bool isintersect = false;
    for(int i =0; i < m.Vertices.size(); ++i){
        if(m.Vertices[i].c == 0){
            vert.push_back(m.Vertices[i]);
        }
    }

    for(int i =0; i < m.Faces.size(); ++i){
        for(int j =0; j < m.Faces[i].Indices.size(); ++j){
            if(m.Vertices[m.Faces[i].Indices[j]].c == 0){
                isintersect = true;
            }
            else {
                isintersect = false;
                break;
            }
        }
        if(isintersect){
            plane.n = Vector{{m.Vertices[m.Faces[i].Indices[0]]},{m.Vertices[m.Faces[i].Indices[1]]}}.cross(Vector{{m.Vertices[m.Faces[i].Indices[1]]},{m.Vertices[m.Faces[i].Indices[2]]}}).normalize();
            plane.p = m.Vertices[m.Faces[i].Indices[0]];
            break;
        }
    }

    if(vert.size() != 0){
        vert = tries(vert, plane);
        //std::cout<<"vert size   "<<vert.size()<<std::endl;
        for(int i = 0; i < vert.size(); ++i){
            //std::cout<<vert[i]<<' '<<' ';
            intersect_face.Indices.push_back(getVertexIndex(vert[i], m));
        }
    }//std::cout<<"vert size   "<<vert.size()<<std::endl;
 



    Mesh m_new;
    std::array<std::vector<Vertex>, 6> arr;
    std::array<Flat,6> p;
    p[0].n = {{1}, {0}, {0}};//это нужно для того чтобы ниже передать в tries
    p[0].p = {{1}, {1}, {1}};

    p[1].n = {{0}, {1}, {0}};
    p[1].p = {{1}, {1}, {1}};

    p[2].n = {{0}, {0}, {1}};
    p[2].p = {{1}, {1}, {1}};

    p[3].n = {{-1}, {0}, {0}};
    p[3].p = {{-1}, {-1}, {-1}};

    p[4].n = {{0}, {-1}, {0}};
    p[4].p = {{-1}, {-1}, {-1}};

    p[5].n = {{0}, {0}, {-1}};
    p[5].p = {{-1}, {-1}, {-1}};

    for(int i =0; i < m.Vertices.size(); ++i){
        if(m.Vertices[i].x == 1){
            arr[0].push_back(m.Vertices[i]);
        }
        if(m.Vertices[i].y == 1){
            arr[1].push_back(m.Vertices[i]);
        }
        if(m.Vertices[i].z == 1){
            arr[2].push_back(m.Vertices[i]);
        }
        if(m.Vertices[i].x == -1){
            arr[3].push_back(m.Vertices[i]);
        }
        if(m.Vertices[i].y == -1){
            arr[4].push_back(m.Vertices[i]);
        }
        if(m.Vertices[i].z == -1){
            arr[5].push_back(m.Vertices[i]);
        }
    }

    for(int i = 0; i < arr.size(); ++i){
        if(arr[i].size() != 0 ){
            Face f_new;

            arr[i] = tries(arr[i], p[i] );
            for(int j =0; j < arr[i].size(); ++j)
                f_new.Indices.push_back(getVertexIndex(arr[i][j], m));
            m_new.Faces.push_back(f_new);
        }
        
    }
    if(intersect_face.Indices.size() >= 3)
        m_new.Faces.push_back(intersect_face);

    m.Faces = m_new.Faces;


}

bool Check(Mesh const &m){
    std::cout<<"v: "<<m.Vertices.size() <<" f:  "<<m.Faces.size()<<"  res:  "<< m.Vertices.size() - m.Faces.size()/2<<std::endl;

    if( m.Vertices.size() - m.Faces.size()/2 == 2)
        return true;
    else
        return false;
}

void Intersect(Mesh &m, Flat const &f){

    Correct(m, f);
    PointClassify(m, f);
    SpecialCases(m, f);

    Mesh res = ResultOfIntersect(m, f);

    Triangulation(res);
    m = res;
    Check(m);

}

