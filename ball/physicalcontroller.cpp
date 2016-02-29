#include "physicalcontroller.h"




/* Detect the ball collision with the ball */
void PhysicalController::findBallBallCollision(Ball *ball1,
                                               Ball *ball2,
                                               double dt,
                                               double x,
                                               GMlib::Array<CollisionObject> &colisionArray)
{
    GMlib::Point<float,3> pointBall1 = ball1->getPos();
    GMlib::Point<float,3> pointBall2 = ball2->getPos();
    GMlib::Vector<float,3> ds1 = ball1->getds();
    GMlib::Vector<float,3> ds2 = ball2->getds();

    double t1,t2;
    t1 = ball1->getX();
    t2 = ball2->getX();

    GMlib::Vector<float,3> ds = (1+t1)*ds1 - (1+t2)*ds2;
    GMlib::Vector<float,3> p1 = pointBall1 - t1*ds1 ;
    GMlib::Vector<float,3> p2 = pointBall2 - t2*ds2;
    GMlib::Vector<float,3> dp = p1-p2 ;

    float r = ball1->getRadius() + ball2->getRadius();
    double a = ds*ds;
    double b = 2*ds*dp;
    double c = dp*dp -r*r;
    double delta = b*b-4*a*c;

    //check the ball and ball if intersected
    if((pointBall1-pointBall2).getLength() < ball1->getRadius() + ball2->getRadius())
    {   double s= 0.51*(dp.getLength()-ball1->getRadius()-ball2->getRadius() );

        ball1->translate(s*dp);
        ball2->translate(-s*dp);

        dp *= 1+(2*s);
        b = 2*(dp*ds);
        c = (dp*dp) - (r*r);

    }


    if(delta < 0)
    {
    }
    else if(delta == 0)
    {
    }
    else
    {
        x = (-b-sqrt(delta))/(2*a);
        if(x>0 && x<=1)
        {
            colisionArray.insertAlways(CollisionObject(ball1,ball2,x),true);
        }
    }
}

/* Detect the ball collision with the wall */
void PhysicalController::findBallWallCollision(Ball *ball,
                                               GMlib::PPlane<float> *wall,
                                               double dt,
                                               double x,
                                               GMlib::Array<CollisionObject> &collisionArray)
{
    GMlib::Point<float,3> pointBall = ball->getPos();
    wall->getClosestPoint(this->getPos(),_u,_v);
    GMlib::DMatrix<GMlib::Vector<float,3>> mat=wall->evaluate(_u,_v,1,1);
    float radious = ball->getRadius();

    GMlib::Vector<float,3> d = mat[0][0]-pointBall ;
    GMlib::Vector<float,3> n = wall->getNormal();
    GMlib::Vector<float,3> ds = ball->getds();

    double DN=d*n;

    if(DN+radious>0.0)   //check the ball and wall if intersected
    {
        ball->translate(2.0*(d*n+radious)*wall->getNormal() );
        DN-=2.0*(DN+radious);
    }

    x = (radious+(DN))/(ds*n);
    if(x >0 && x<=1)
    {
        collisionArray.insertAlways(CollisionObject(ball,wall,x));
    }

}



/* set the new velocity of the ball after collide with the wall*/

void PhysicalController::collisionBallWall(Ball *ball, GMlib::PPlane<float> *plane, double dt)
{

    GMlib::Vector<float,3> _ve = ball->getVelocity();
    _ve -= 2*(plane->getNormal()*_ve)*plane->getNormal();  //球与墙相撞后的速度公式

    if(ball->getVelocity().getLength() > 0.1)
    {
    ball->setVelocity(_ve);

    ball->updateStep(dt);
    }
    else
    {
        ball->setVelocity(GMlib::Vector<float,3>(0,0,0));
    }

}


/*  Calculate the velocity of the ball after colliding with another ball */

void PhysicalController::collisionBallBall(Ball *ball1, Ball *ball2, double dt)
{
    GMlib::Vector<float,3> newvelocity1 = ball1->getVelocity();
    GMlib::UnitVector<float,3> d = ball2->getPos() - ball1->getPos();
    GMlib::Vector<float,3> v1 = (newvelocity1*d/(d*d))*d;
    GMlib::Vector<float,3> v1n = newvelocity1 - v1;

    GMlib::Vector<float,3> newvelocity2 = ball2->getVelocity();
    GMlib::Vector<float,3> v2 = (newvelocity2*d/(d*d))*d;
    GMlib::Vector<float,3> v2n = newvelocity2 - v2;

    double M1, M2, M3;
    M1 = (ball1->getMass() - ball2->getMass()) / (ball1->getMass() + ball2->getMass());
    M2 = (2*ball2->getMass()/(ball1->getMass() + ball2->getMass()));
    M3 = (2*ball1->getMass()/(ball1->getMass() + ball2->getMass()));
    GMlib::Vector<float,3> v1after = M1*v1 + M2*v2;
    GMlib::Vector<float,3> v2after = -M1*v2 + M3*v1;

    GMlib::Vector<float,3> newvelocityball1 = v1after + v1n;
    GMlib::Vector<float,3> newvelocityball2 = v2after + v2n;

    ball1->setVelocity(newvelocityball1*0.95);
    ball2->setVelocity(newvelocityball2*0.95);

    ball1->updateStep(dt);
    ball2->updateStep(dt);
}

void PhysicalController::localSimulate(double dt)
{
    updateStep(dt);
    coll_dectect(dt);
}

void PhysicalController::updateStep(double dt)
{
    for(int i=0;i<_ballArray.getSize();i++)
        _ballArray[i]->updateStep(dt);
}

void PhysicalController::coll_dectect(double dt)
{
    {
        for (int i=0;i<_ballArray.size();i++)
        {
            for(int j=i+1;j<_ballArray.size();j++)
            {
                findBallBallCollision(_ballArray[i],_ballArray[j],dt,0,_collisionArray);
            }
        }
    }

    for (int i=0;i<_ballArray.size();i++)
    {
        for(int j=0;j<_planeArray.size();j++)
        {
            findBallWallCollision(_ballArray[i],_planeArray[j],dt,0,_collisionArray);
        }
    }

    while(_collisionArray.getSize()>0)
    {
        _collisionArray.sort();
        _collisionArray.makeUnique();
        CollisionObject collisionObject = _collisionArray[0];
        _collisionArray.removeIndex(0);

        if(collisionObject.isBW())
        {
            collisionBallWall(collisionObject.getBall(0),
                              collisionObject.getPlane(),
                              (1-collisionObject.getX())*dt);
        }
        else
        {
            collisionBallBall(collisionObject.getBall(0),
                              collisionObject.getBall(1),
                              (1-collisionObject.getX())*dt);
        }
    }
}


