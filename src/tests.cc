#include "tests.h"





void Test(Mesh const &m_in, double precise){
    Mesh m{m_in};

    Flat plane2;

    plane2.p = {{1}, {1}, {0.99994}}; // empty
    plane2.n = Vector{{0}, {0}, {-1}}.normalize();
    Intersect(m, plane2, precise);
    if(!Check( m ))
        std::cout<<"OK   "<<  Check(m)<<std::endl;
    else
        std::cout<<" ne OK"<<std::endl;
    m = m_in;

    

    plane2.p = {{1}, {1}, {0.99993}}; // not empty
    Intersect(m, plane2, precise);
    if(Check( m ))
        std::cout<<"OK   "<<  Check(m)<<std::endl;
    else
        std::cout<<" ne OK"<<std::endl;
    m = m_in;


    plane2.p = {{1}, {1}, {0.9999}};// empty
    plane2.n = Vector{{-1}, {-1}, {-1}}.normalize();
    Intersect(m, plane2, precise);
    if(!Check( m ))
        std::cout<<"OK   "<<  Check(m)<<std::endl;
    else
        std::cout<<" ne OK"<<std::endl;
    m = m_in;

    plane2.p = {{1}, {1}, {0.9998}};//not empty
    Intersect(m, plane2, precise);
    if(Check( m ))
        std::cout<<"OK   "<<  Check(m)<<std::endl;
    else
        std::cout<<" ne OK"<<std::endl;
    m = m_in;

    plane2.p = {{0}, {0}, {1}}; // empty
    plane2.n = Vector{{0}, {0.00006}, {-1}}.normalize();
    Intersect(m, plane2, precise);
    if(!Check( m))
        std::cout<<"OK   "<<  Check(m)<<std::endl;
    else
        std::cout<<" ne OK"<<std::endl;
    m = m_in;
                                 // not empty
    plane2.n = Vector{{0}, {0.00007}, {-1}}.normalize();
    Intersect(m, plane2, precise);
    if(Check( m))
        std::cout<<"OK   "<<  Check(m)<<std::endl;
    else
        std::cout<<" ne OK"<<std::endl;

}


void Icosahedron( Mesh &m, double precise){
    double a = 0.525731112;
    double b = 0.850650808;
    Vertex v0 ={ -a , 0 , b};
    Vertex v1 ={ a , 0 , b};
    Vertex v2 ={ -a , 0 , -b};
    Vertex v3 ={ a , 0 , -b};
    Vertex v4 ={ 0, b, a};
    Vertex v5 ={ 0, b, -a};
    Vertex v6 ={ 0, -b, a};
    Vertex v7 ={ 0, -b, -a};
    Vertex v8 ={ b, a, 0};
    Vertex v9 ={ -b, a, 0};
    Vertex v10 ={ b, -a, 0};
    Vertex v11 ={ -b, -a, 0};

    Intersect(m, FlatByPoints(v0, v4, v1), precise);
    Intersect(m, FlatByPoints(v0, v9, v4), precise);
    Intersect(m, FlatByPoints(v9, v5, v4), precise);
    Intersect(m, FlatByPoints(v4, v5, v8), precise);
    Intersect(m, FlatByPoints(v4, v8, v1), precise);

    Intersect(m, FlatByPoints(v8, v10, v1), precise);
    Intersect(m, FlatByPoints(v8, v3, v10), precise);
    Intersect(m, FlatByPoints(v5, v3, v8), precise);
    Intersect(m, FlatByPoints(v5, v2, v3), precise);
    Intersect(m, FlatByPoints(v2, v7, v3), precise);

    Intersect(m, FlatByPoints(v7, v10, v3), precise);
    Intersect(m, FlatByPoints(v7, v6, v10), precise);
    Intersect(m, FlatByPoints(v7, v11, v6), precise);
    Intersect(m, FlatByPoints(v11, v0, v6), precise);
    Intersect(m, FlatByPoints(v0, v1, v6), precise);

    Intersect(m, FlatByPoints(v6, v1, v10), precise);
    Intersect(m, FlatByPoints(v9, v0, v11), precise);
    Intersect(m, FlatByPoints(v9, v11, v2), precise);
    Intersect(m, FlatByPoints(v9, v2, v5), precise);
    Intersect(m, FlatByPoints(v7, v2, v11), precise);

 //cut

    v0 ={ -a *1.5 , 0 , b*1.5};
    v1 ={ a*1.5 , 0 , b*1.5};
    v2 ={ -a*1.5 , 0 , -b*1.5};
    v3 ={ a*1.5 , 0 , -b*1.5};
    v4 ={ 0, b*1.5, a*1.5};
    v7 ={ 0, -b*1.5, -a*1.5};
    v8 ={ b*1.5, a*1.5, 0};

    Intersect(m, FlatByPoints(v7, v2, v0), precise);
    Intersect(m, FlatByPoints(v7, v1, v3), precise);
    Intersect(m, FlatByPoints(v0, v8, v1), precise);



    // Intersect(m, FlatByPoints(v4, v2, v10), precise);
    // Intersect(m, FlatByPoints(v6, v2, v1), precise);

    // Intersect(m, FlatByPoints(v5, v11, v10), precise);

}

void Cuboctahedron(Mesh &m, double precise){
    
    Flat plane;

    plane.p = {{0.5}, {0.5}, {1}};
    plane.n = Vector{{1}, {1}, {1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{-0.5}, {-0.5}, {1}};
    plane.n = Vector{{-1}, {-1}, {1}}.normalize();
    Intersect(m, plane, precise);
    
    plane.p = {{-0.5}, {0.5}, {-1}};
    plane.n = Vector{{-1}, {1}, {-1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{0.5}, {-0.5}, {-1}};
    plane.n = Vector{{1}, {-1}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{-0.5}, {0.5}, {1}};
    plane.n = Vector{{-1}, {1}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{0.5}, {-0.5}, {1}};
    plane.n = Vector{{1}, {-1}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{0.5}, {0.5}, {-1}};
    plane.n = Vector{{1}, {1}, {-1}}.normalize();
    Intersect(m, plane, precise);



    plane.p = {{-0.5}, {-0.5}, {-1}};
    plane.n = Vector{{-1}, {-1}, {-1}}.normalize();
    Intersect(m, plane, precise);


}

void Rhombicuboctahedron2(Mesh &m, double precise){
    Flat plane;

    plane.p = {{0.3333333333}, {-1}, {1}};
    plane.n = Vector{{1}, {0}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{-1}, {1}, {-0.3333333333}};
    plane.n = Vector{{-1}, {0}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{1}, {0.3333333333}, {1}};
    plane.n = Vector{{0}, {1}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{1}, {-1}, {-0.3333333333}};
    plane.n = Vector{{0}, {-1}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{-0.3333333333}, {1}, {1}};
    plane.n = Vector{{-1}, {0}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{1}, {1}, {-0.3333333333}};
    plane.n = Vector{{1}, {0}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{1}, {-0.3333333333}, {1}};
    plane.n = Vector{{0}, {-1}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{1}, {1}, {-0.3333333333}};
    plane.n = Vector{{0}, {1}, {-1}}.normalize();
    Intersect(m, plane, precise);




    plane.p = {{0.3333333333}, {1}, {1}};
    plane.n = Vector{{1}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{-0.3333333333}, {-1}, {1}};
    plane.n = Vector{{-1}, {-1}, {0}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{0.3333333333}, {-1}, {1}};
    plane.n = Vector{{1}, {-1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{-0.3333333333}, {1}, {1}};
    plane.n = Vector{{-1}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);
}

void Rhombicuboctahedron3(Mesh &m, double precise){
    Flat plane;

    plane.p = {{0.5}, {-1}, {1}};
    plane.n = Vector{{1}, {0}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{-1}, {1}, {-0.5}};
    plane.n = Vector{{-1}, {0}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{1}, {0.5}, {1}};
    plane.n = Vector{{0}, {1}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{1}, {-1}, {-0.5}};
    plane.n = Vector{{0}, {-1}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{-0.5}, {1}, {1}};
    plane.n = Vector{{-1}, {0}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{1}, {1}, {-0.5}};
    plane.n = Vector{{1}, {0}, {-1}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{1}, {-0.5}, {1}};
    plane.n = Vector{{0}, {-1}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{1}, {1}, {-0.5}};
    plane.n = Vector{{0}, {1}, {-1}}.normalize();
    Intersect(m, plane, precise);




    plane.p = {{0.5}, {1}, {1}};
    plane.n = Vector{{1}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{-0.5}, {-1}, {1}};
    plane.n = Vector{{-1}, {-1}, {0}}.normalize();
    Intersect(m, plane, precise);


    plane.p = {{0.5}, {-1}, {1}};
    plane.n = Vector{{1}, {-1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{-0.5}, {1}, {1}};
    plane.n = Vector{{-1}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);
}

void Rhombicuboctahedron(Mesh &m, double precise){
    Cuboctahedron( m , precise); 
    Rhombicuboctahedron3( m , precise);
}

void Pyramid(Mesh &m, double precise){
    Flat plane;

    plane.p = {{0}, {1}, {0}};

    plane.n = Vector{{2}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{-2}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{0}, {1}, {2}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{0}, {1}, {-2}}.normalize();
    Intersect(m, plane, precise);

}

void Octahedron(Mesh &m, double precise){
    Flat plane;

    plane.p = {{0}, {-1}, {0}};

    plane.n = Vector{{1.144213}, {-1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{-1.144213}, {-1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{0}, {-1}, {1.144213}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{0}, {-1}, {-1.144213}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{0}, {1}, {0}};

    plane.n = Vector{{1.144213}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{-1.144213}, {1}, {0}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{0}, {1}, {1.144213}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{0}, {1}, {-1.144213}}.normalize();
    Intersect(m, plane, precise);

}

void Tetrahedron(Mesh &m, double precise){
    Flat plane;

    plane.p = {{0}, {0}, {1}};

    plane.n = Vector{{1.144213}, {0}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{-1}, {1.744213}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.n = Vector{{-1}, {-1.744213}, {1}}.normalize();
    Intersect(m, plane, precise);

    plane.p = {{0}, {0}, {0.5}};
    plane.n = Vector{{0}, {0}, {-1}}.normalize();
    Intersect(m, plane, precise);
}


