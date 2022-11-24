#include "structures.h"
#include "obj_parser.h"


float ScalarProduct(Vector &r1, Vector &r2){
    return (r1.x * r2.x + r1.y * r2.y + r1.z * r2.z);
};

Vector VectorProduct(Vector r1, Vector r2){
    Vector r;
    r.x = (r1.y * r2.z) - (r1.z * r2.y);
    r.y = -(r1.x * r2.z) + (r1.z * r2.x);
    r.z = (r1.x * r2.y) - (r1.y * r2.x);
    return r;
};

int PointInFlat (Vertex &p, Flat &f){//подставление координат точки в ур.плоскости
    return f.A * p.x + f.B * p.y + f.C * p.z + f.D;
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

    std::array arr[intersect.size()-1];//массив - количество появления вершин на второй позиции
    for(int i=0; i<intersect.size() - 1; ++i)
        arr[i] = 0;
    for( int i = 1; i < intersect.size() - 2; ++i ){
        for(int j = i + 1; j < intersect.size(); ++j ){

            Vector try = VectorProduct( Vector(intersect[0], intersect[i]), Vector(intersect[i], intersect[j]));
            try.normalize();
            if(try.x == norm.x && try.y == norm.y && try.z == norm.z ){
                arr[i-1] ++;
            }
            else if(try.x == -norm.x && try.y == -norm.y && try.z == -norm.z ){
                arr[j-1] ++;
            }  
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

Vertex Segment_Flat_Intersection(Segment s, Flat f){//точка пересечения плоскости и отрезка
    if(PointInFlat(s.A,f) * PointInFlat(s.B,f) < 0){
        Vertex p = operator-(s.A, s.B);
        float t = - ( PointInFlat(s.B, f)) / (PointInFlat(p, f) - f.D);
        return operator+( operator* (s.A,t), operator* (s.B,1-t));
    }
}

int getVertexIndex(Vertex v, Mesh m){
    int m = m.Vertices.size();
    for(int i =0; i < m.Vertices.size(); ++i){
        m = m-1;
        if(operator==(m.Vertices[i], v ))
            return i;
    }
    if (m == 0){
        return;
    }
}

//не работает, если плоскость пересекает вершины корректно

void ResultOfIntersect( Mesh &m, Flat &f){
    int index = m.Vertices.size();//индексы новых вершин
    int zero = OnPoints(m, f);//количество вершин объекта, лежащих на плоскости
    int sum = InPoints(m, f) - OutPoints(m, f);//cумма codes. нужно для определения крайних/вырожденных случаев

    Face new_face;
    m.Faces.push_back(new_face);
    
    if((sum == - (m.Vertices.size() - zero)) || (zero == 4 && sum == -4 )){//случаи, где результат - пустое множество
        std::cerr << "Empty intersect " << std::endl;
        return;
    }
    if((sum == (m.Vertices.size() - zero)) || (zero == 4 && sum == 4 )){// случаи, когда результат - исходный объект
        return m;
    }

    std::vector<Vertex> intersect; //points of intersect
    for(int i = 0; i < m.Faces.size(); ++i){//для каждой грани

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){//обходим каждую вершину грани

            if( m.Vertices[m.Faces[i].Indices[j]].c != m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]].c ){ //if adjacent vertices with different codes
                Segment s;
                s.A = m.Vertices[m.Faces[i].Indices[j]];
                s.B = m.Vertices[m.Faces[i].Indices[(j+1) % m.Faces[i].Indices.size()]];
                Vertex r = Segment_Flat_Intersection(s, f);
                r.c = 0;
                intersect.push_back(r);
                index ++;

                //вставка индекса новой вершины в конец списвка грани
                if(j == m.Faces[i].Indices.size() - 1){
                    m.Faces[i].Indices.push_back(index);
                    //еще вставить в MESH
                }
                else{
                //вставка индекса новой вершины в список на место j+1 (после j) в списке грани
                    m.Faces[i].Indices.push_back( m.Faces[i].Indices[m.Faces[i].Indices.size() - 1] );//копируем последний элемент в конец
                    for(int s = m.Faces[i].Indices.size() - 2; s > j; s--){
                        m.Faces[i].Indices[s] = m.Faces[i].Indices[s-1];
                    }
                    m.Faces[i].Indices[j+1] = index;
                    //еще вставить в MESH
                }
            }
        }

        for(int j =0; j < m.Faces[i].Indices.size(); ++j){ //удаление из списка вершины//////////////////////////////
            if( m.Vertices[m.Faces[i].Indices[j]].c == -1){
                for(int s = j; s < m.Faces[i].Indices.size() - 1 ; ++s){
                    m.Faces[i].Indices[s] = m.Faces[i].Indices[s+1];
                }
                m.Faces[i].Indices.pop_back;//удалить с конца?
                //еще удалить в MESH
            }
        }

    }
    
    intersect = tries(intersect, f);//вектор вершин в нужном порядке для новой грани
    for(int i = 0; i<intersect.size(); ++i){
        new_face.Indices.push_back(intersect[i]);
    }
    
    m.Faces.push_back(new_face);





    //запись в файл новой грани и новых вершин
   


}