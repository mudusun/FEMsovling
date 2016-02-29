#include "collisionobject.h"


 /*construction of the collisionobject class with ball and plane parameter */

CollisionObject::CollisionObject(Ball *ball, GMlib::PPlane<float> *plane, double x)
{
  this->ball[0] = ball;
  this->wall = plane;
  this->x = x;
  this->bw = true;
}


 /* construction of the CollisionObject class with ball and ball parameter */

CollisionObject::CollisionObject(Ball *ball1, Ball *ball2, double x)
{
  this->ball[0] = ball1;
  this->ball[1] = ball2;
  this->x = x;
  this->bw = false;
}


bool CollisionObject::operator == (const CollisionObject& other) const
{
  if(ball[0] == other.ball[0])
    return true;
  if(!other.bw && ball[0] == other.ball[1])
    return true;
  if(ball[1] == other.ball[0])
    return true;
  if(!other.bw && ball[1] == other.ball[1])
    return true;
  return false;
}


/*the construction function with no parameter*/

CollisionObject::CollisionObject()
{

}

/*this function is the override function in the operator */
bool CollisionObject::operator < (const CollisionObject& other) const
{
  return x < other.x;
}

/*get the plane*/
GMlib::PPlane<float>* CollisionObject::getPlane() const
{
  return wall;
}

/* get the ball*/
Ball* CollisionObject::getBall(int i)
{
  return ball[i];
}

/*get the x value */
double CollisionObject::getX()
{
  return x;
}

/* return the value ball and wall collision or not */
bool CollisionObject::isBW()
{
  return bw;
}
