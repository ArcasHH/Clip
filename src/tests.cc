#include "tests.h"


void Cuboctahedron(Mesh &m){
    
    Flat plane;

    plane.p = {{0.5}, {0.5}, {1}};
    plane.n = Vector{{1}, {1}, {1}}.normalize();
    Intersect(m, plane);


    plane.p = {{-0.5}, {-0.5}, {1}};
    plane.n = Vector{{-1}, {-1}, {1}}.normalize();
    Intersect(m, plane);
    
    plane.p = {{-0.5}, {0.5}, {-1}};
    plane.n = Vector{{-1}, {1}, {-1}}.normalize();
    Intersect(m, plane);

    plane.p = {{0.5}, {-0.5}, {-1}};
    plane.n = Vector{{1}, {-1}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{-0.5}, {0.5}, {1}};
    plane.n = Vector{{-1}, {1}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{0.5}, {-0.5}, {1}};
    plane.n = Vector{{1}, {-1}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{0.5}, {0.5}, {-1}};
    plane.n = Vector{{1}, {1}, {-1}}.normalize();
    Intersect(m, plane);



    plane.p = {{-0.5}, {-0.5}, {-1}};
    plane.n = Vector{{-1}, {-1}, {-1}}.normalize();
    Intersect(m, plane);


}

void Rhombicuboctahedron2(Mesh &m){
    Flat plane;

    plane.p = {{0.3333333333}, {-1}, {1}};
    plane.n = Vector{{1}, {0}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{-1}, {1}, {-0.3333333333}};
    plane.n = Vector{{-1}, {0}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{1}, {0.3333333333}, {1}};
    plane.n = Vector{{0}, {1}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{1}, {-1}, {-0.3333333333}};
    plane.n = Vector{{0}, {-1}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{-0.3333333333}, {1}, {1}};
    plane.n = Vector{{-1}, {0}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{1}, {1}, {-0.3333333333}};
    plane.n = Vector{{1}, {0}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{1}, {-0.3333333333}, {1}};
    plane.n = Vector{{0}, {-1}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{1}, {1}, {-0.3333333333}};
    plane.n = Vector{{0}, {1}, {-1}}.normalize();
    Intersect(m, plane);




    plane.p = {{0.3333333333}, {1}, {1}};
    plane.n = Vector{{1}, {1}, {0}}.normalize();
    Intersect(m, plane);

    plane.p = {{-0.3333333333}, {-1}, {1}};
    plane.n = Vector{{-1}, {-1}, {0}}.normalize();
    Intersect(m, plane);


    plane.p = {{0.3333333333}, {-1}, {1}};
    plane.n = Vector{{1}, {-1}, {0}}.normalize();
    Intersect(m, plane);

    plane.p = {{-0.3333333333}, {1}, {1}};
    plane.n = Vector{{-1}, {1}, {0}}.normalize();
    Intersect(m, plane);
}

void Rhombicuboctahedron3(Mesh &m){
    Flat plane;

    plane.p = {{0.5}, {-1}, {1}};
    plane.n = Vector{{1}, {0}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{-1}, {1}, {-0.5}};
    plane.n = Vector{{-1}, {0}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{1}, {0.5}, {1}};
    plane.n = Vector{{0}, {1}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{1}, {-1}, {-0.5}};
    plane.n = Vector{{0}, {-1}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{-0.5}, {1}, {1}};
    plane.n = Vector{{-1}, {0}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{1}, {1}, {-0.5}};
    plane.n = Vector{{1}, {0}, {-1}}.normalize();
    Intersect(m, plane);


    plane.p = {{1}, {-0.5}, {1}};
    plane.n = Vector{{0}, {-1}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{1}, {1}, {-0.5}};
    plane.n = Vector{{0}, {1}, {-1}}.normalize();
    Intersect(m, plane);




    plane.p = {{0.5}, {1}, {1}};
    plane.n = Vector{{1}, {1}, {0}}.normalize();
    Intersect(m, plane);

    plane.p = {{-0.5}, {-1}, {1}};
    plane.n = Vector{{-1}, {-1}, {0}}.normalize();
    Intersect(m, plane);


    plane.p = {{0.5}, {-1}, {1}};
    plane.n = Vector{{1}, {-1}, {0}}.normalize();
    Intersect(m, plane);

    plane.p = {{-0.5}, {1}, {1}};
    plane.n = Vector{{-1}, {1}, {0}}.normalize();
    Intersect(m, plane);
}


void Pyramid(Mesh &m){
    Flat plane;

    plane.p = {{0}, {1}, {0}};

    plane.n = Vector{{2}, {1}, {0}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{-2}, {1}, {0}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{0}, {1}, {2}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{0}, {1}, {-2}}.normalize();
    Intersect(m, plane);

}

void Octahedron(Mesh &m){
    Flat plane;

    plane.p = {{0}, {-1}, {0}};

    plane.n = Vector{{1.144213}, {-1}, {0}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{-1.144213}, {-1}, {0}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{0}, {-1}, {1.144213}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{0}, {-1}, {-1.144213}}.normalize();
    Intersect(m, plane);

    plane.p = {{0}, {1}, {0}};

    plane.n = Vector{{1.144213}, {1}, {0}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{-1.144213}, {1}, {0}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{0}, {1}, {1.144213}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{0}, {1}, {-1.144213}}.normalize();
    Intersect(m, plane);

}

void Tetrahedron(Mesh &m){
     Flat plane;

    plane.p = {{0}, {0}, {1}};

    plane.n = Vector{{1.144213}, {0}, {1}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{-1}, {1.144213}, {1}}.normalize();
    Intersect(m, plane);

    plane.n = Vector{{-1}, {-1.144213}, {1}}.normalize();
    Intersect(m, plane);

    plane.p = {{0}, {0}, {0.5}};
    plane.n = Vector{{0}, {0}, {-1}}.normalize();
    Intersect(m, plane);
}