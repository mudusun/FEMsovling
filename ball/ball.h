/*
 * This class is define the ball, this class include
 * _mass: the mass of the ball
 * _x: the value to store the time
 * _Surf: the plane object
 * _p: the position of the ball
 * _q: the point of the plane
 * _velocity: the velocity of the ball
 * _u: the u coordinate
 * _v: the v coordinate
 * _normal: the normal of the wall
 * _ds: the ds vector of the ball
*/

#ifndef BALL_H
#define BALL_H

#include <parametrics/gmpsphere>
#include "plane.h"

class Ball : public GMlib::PSphere<float> {

private:
  float _mass;
  double _x;
  Plane<float>* _Surf;
  GMlib::Point<float,3> _p;
  GMlib::Point<float,3> _q;
  GMlib::Vector<float,3> _velocity;
  GMlib::UnitVector<float,3> _normal;
  GMlib::Vector<float,3> _ds;
  float _u, _v;
  bool _go;

public:
  using PSphere::PSphere;
  Ball(GMlib::Vector<float,3> velocity, float mass, float radious, Plane<float>* surf):GMlib::PSphere<float>(radious)
  {
    this->_velocity = velocity;
    this->_Surf = surf;
    this->_mass = mass;
    _x = 0;
    _go = false;
  }
  ~Ball();

  void updateStep(double dt);
  void setVelocity(GMlib::Vector<float,3> newVelocity);

  void setX(double x);

  GMlib::Vector<float,3> getds();
  GMlib::Vector<float,3> getSurfNormal();
  GMlib::Vector<float,3> getVelocity();

  float  getMass();
  double getX();


  void moveUp();
  void moveDown();
  void moveRight();
  void moveLeft();

protected:
  void localSimulate(double dt) override ;

};

#endif


