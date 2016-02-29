#ifndef COLLISIONOBJECT_H
#define COLLISIONOBJECT_H

#include "ball.h"
#include "plane.h"
#include <parametrics/surfaces/gmpplane.h>
class CollisionObject{

private:
  Ball *ball[2];               //define the ball array
  GMlib::PPlane<float> *wall;  //define the wall
  double x;                    //time
  bool bw;                     //check the ball and wall are collided or not

public:
  CollisionObject();
  CollisionObject(Ball *ball1, Ball *ball2, double x);
  CollisionObject(Ball *ball, GMlib::PPlane<float> *plane, double x);

  bool operator < (const CollisionObject& other) const;
  bool operator == (const CollisionObject& other) const;

  GMlib::PPlane<float>* getPlane() const;
  Ball* getBall(int i);
  double getX();

  bool isBW();

};

#endif
