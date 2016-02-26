#ifndef ELEMENT_H
#define ELEMENT_H
#include <trianglesystem/gmtrianglesystem>
#include <core/gmarraylx>
#include <math.h>
#include "element.h"

class element
{

GMlib::TSVertex<float> getVertex();

public:
    element(){}
    element(GMlib::TSVertex<float> &vertex)
    {
        node=&vertex;
    }


    ~element(){}

void findElement();

GMlib::Array<GMlib::TSEdge<float>*> getEdges()
{
    return node->getEdges();
}

GMlib::Array<GMlib::TSTriangle<float>*> getTriangles()
{
    return node->getTriangles();
}




private:
GMlib::TSVertex<float>* node;


};

#endif // ELEMENT_H

