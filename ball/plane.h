#ifndef PLANE_H
#define PLANE_H


#include <QTimerEvent>
#include <QRectF>
#include <QMouseEvent>
#include <QDebug>
#include <stdexcept>
#include <thread>
#include <mutex>
#include <parametrics/gmpplane>
#include <parametrics/gmpbeziersurf>

template <typename T>
class Plane: public GMlib::PSurf<T,3>
{
    GM_SCENEOBJECT(Plane)

    public:

        Plane( const GMlib::DMatrix<GMlib::Vector<T,3> >& m);
        Plane( const Plane<T>& copy );
    virtual ~Plane();

protected:
    GMlib::DMatrix<GMlib::Vector<T,3> >	_m;

    void                      eval(T u, T v, int d1, int d2, bool lu = true, bool lv = true );
    T                         getEndPU();
    T                         getEndPV();
    T                         getStartPU();
    T                         getStartPV();
};


template <typename T>
inline
Plane<T>::Plane( const GMlib::DMatrix<GMlib::Vector<T,3> >& m )
{
    _m = m;
    this->_dm = GMlib::GM_DERIVATION_EXPLICIT;
}


template <typename T>
inline
Plane<T>::Plane( const Plane<T>& copy ):GMlib::PSurf<T,3>(copy)
{
    _m   = copy._m;
}


template <typename T>
Plane<T>::~Plane() {}

template <typename T>
void Plane<T>::eval(T u, T v, int d1, int d2, bool lu, bool lv )
{

    this->_p.setDim( d1+1, d2+1 );

    GMlib::DVector<T> u1(4);
    u1[0]=(1-u)*(1-u)*(1-u);
    u1[1]=3*u*(1-u)*(1-u);
    u1[2]=3*u*u*(1-u);
    u1[3]=u*u*u;
    GMlib::DVector<T> v1(4);
    v1[0]=(1-v)*(1-v)*(1-v);
    v1[1]=3*v*(1-v)*(1-v);
    v1[2]=3*v*v*(1-v);
    v1[3]=v*v*v;

    this->_p[0][0] = u1*(_m^v1);
    if( this->_dm == GMlib::GM_DERIVATION_EXPLICIT )
    {
        // 1st
        GMlib::DVector<T>uu1(4);
        uu1[0]=(-3*(1-u)*(1-u));
        uu1[1]=3*(1-u)*(1-3*u);
        uu1[2]=6*u-9*u*u;
        uu1[3]=3*u*u;
        GMlib::DVector<T>vv1(4);
        vv1[0]=(-3*(1-v)*(1-v));
        vv1[1]=3*(1-v)*(1-3*v);
        vv1[2]=6*v-9*v*v;
        vv1[3]=3*v*v;

        if(d1)            this->_p[1][0] = uu1*(_m^v1); // S_u
        if(d2)            this->_p[0][1] = u1*(_m^vv1); // S_v
        if(d1>1 && d2>1)  this->_p[1][1] = uu1*(_m^vv1); // S_uv


        // 2nd
        if(d1>1)          this->_p[2][0] = GMlib::Vector<T,3>(T(0)); // S_uu
        if(d2>1)          this->_p[0][2] = GMlib::Vector<T,3>(T(0)); // S_vv
        if(d1>1 && d2)    this->_p[2][1] = GMlib::Vector<T,3>(T(0)); // S_uuv
        if(d1   && d2>1)  this->_p[1][2] = GMlib::Vector<T,3>(T(0)); // S_uvv
        if(d1>1 && d2>1)  this->_p[2][2] = GMlib::Vector<T,3>(T(0)); // S_uuvv

        // 3rd
        if(d1>2)          this->_p[3][0] = GMlib::Vector<T,3>(T(0)); // S_uuu
        if(d2>2)          this->_p[0][3] = GMlib::Vector<T,3>(T(0)); // S_vvv
        if(d1>2 && d2)    this->_p[3][1] = GMlib::Vector<T,3>(T(0)); // S_uuuv
        if(d1   && d2>2)  this->_p[1][3] = GMlib::Vector<T,3>(T(0)); // S_uvvv
        if(d1>2 && d2>1)  this->_p[3][2] = GMlib::Vector<T,3>(T(0)); // S_uuuvv
        if(d1>1 && d2>2)  this->_p[2][3] = GMlib::Vector<T,3>(T(0)); // S_uuvvv
        if(d1>2 && d2>2)  this->_p[3][3] = GMlib::Vector<T,3>(T(0)); // S_uuuvvv

        // 4th
        if(d1>3)          this->_p[4][0] = GMlib::Vector<T,3>(T(0)); // S_uuuu
        if(d2>3)          this->_p[0][4] = GMlib::Vector<T,3>(T(0)); // S_vvvv
        if(d1>3 && d2)    this->_p[4][1] = GMlib::Vector<T,3>(T(0)); // S_uuuuv
        if(d1   && d2>3)  this->_p[1][4] = GMlib::Vector<T,3>(T(0)); // S_uvvvv
        if(d1>3 && d2>1)  this->_p[4][2] = GMlib::Vector<T,3>(T(0)); // S_uuuuvv
        if(d1>1 && d2>3)  this->_p[2][4] = GMlib::Vector<T,3>(T(0)); // S_uuvvvv
        if(d1>3 && d2>2)  this->_p[4][3] = GMlib::Vector<T,3>(T(0)); // S_uuuuvvv
        if(d1>2 && d2>3)  this->_p[3][4] = GMlib::Vector<T,3>(T(0)); // S_uuuvvvv
        if(d1>3 && d2>3)  this->_p[4][4] = GMlib::Vector<T,3>(T(0)); // S_uuuuvvvv
    }
}


template <typename T>
inline
T Plane<T>::getEndPU()
{
    return T(1);
}

template <typename T>
inline
T Plane<T>::getEndPV()
{
    return T(1);
}


template <typename T>
inline
T Plane<T>::getStartPU()
{
    return T(0);
}


template <typename T>
inline
T Plane<T>::getStartPV()
{
    return T(0);
}

#endif
