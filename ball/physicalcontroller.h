/*
 this is the most importaint class
 _surf: define the plane
 _ballArray: define the array which containt the ball
 _planeArray: define the array which containt the wall
 _collisionArray: store the data of the collision
 */

#ifndef PHYSICALCONTROLLER_H
#define PHYSICALCONTROLLER_H

#include "ball.h"
#include "plane.h"
#include "collisionobject.h"
#include <parametrics/surfaces/gmpplane.h>
#include <core/gmarray>
#include <QTimerEvent>
#include <QRectF>
#include <QMouseEvent>
#include <QDebug>
#include <stdexcept>
#include <thread>
#include <mutex>

class PhysicalController : public GMlib::PSphere<float>
{
    GM_SCENEOBJECT(PhysicalController)

    private:
    Plane<float>* _surf;
//    GMlib::PBezierSurf<float>* _surf;
    GMlib::Array<Ball*> _ballArray;
    GMlib::Array<GMlib::PPlane<float>*> _planeArray;
    GMlib::Array<CollisionObject> _collisionArray;

public:


    PhysicalController(Plane<float>* plane)
//    PhysicalController(GMlib::PBezierSurf<float>* plane)
    {
        this->_surf=plane;
        this->toggleDefaultVisualizer();
        this->replot(20,20,1,1);
        this->setVisible(false);
        this->insert(_surf);
    }

    void addball(Ball* ball)
    {
        this->insert(ball);
        _ballArray+=ball;
    }

    void addplane(GMlib::PPlane<float>* plane)
    {
        this->insert(plane);
        _planeArray+=plane;

    }
    void insertObject();

protected:

    void collisionBallWall(Ball *ball, GMlib::PPlane<float> *plane, double dt);
    void collisionBallBall(Ball *ball1, Ball *ball2, double dt);

    void localSimulate(double dt);

    void updateStep(double dt);
    void coll_dectect(double dt);

    void findBallBallCollision(Ball *ball1,
                               Ball *ball2,
                               double dt,
                               double x,
                               GMlib::Array<CollisionObject> &colisionArray);
    void findBallWallCollision(Ball *ball,
                               GMlib::PPlane<float> *wall,
                               double dt,
                               double x,
                               GMlib::Array<CollisionObject>& collisionArray);

};


#endif
