#ifndef MESH_H
#define MESH_H

#include <trianglesystem/gmtrianglesystem>
#include <core/gmarraylx>
#include <math.h>

#include "element.h"


using namespace std;

class mesh : public GMlib::TriangleFacets<float> {

public:
    mesh()
    {
        force=0;
        swp=true;
    }
    ~mesh(){}

    void random(float rad, int triangle);
    GMlib::TSEdge<float>* isNeibourPoints(GMlib::TSVertex<float>* n0,GMlib::TSVertex<float>* n1);
    void findElement();
    float gen_float(float a, float b);
    void regular(int circle, double firstRad,int firstNode);
    GMlib::Vector<GMlib::Vector<float,2>,3> findVector( GMlib::TSEdge<float>* commonEdge);
    GMlib::Vector<GMlib::Vector<float,2>,3> findVectord(GMlib::TSVertex<float>* point,
                                                        GMlib::TSTriangle<float>* triangle);

private:

    GMlib::ArrayLX<GMlib::TSVertex<float>*> nodes;

    GMlib::DMatrix<float> A;
    GMlib::DVector<float> b;
    GMlib::DVector<float> x;
    GMlib::DMatrix<float> Ainverted;
    bool checkP(GMlib::TSVertex<float>* p, float s, float radius);
    float force;
    bool swp;

protected:
    void localSimulate(double dt);


};





#endif // MESH_H
