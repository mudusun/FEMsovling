#include "mesh.h"
#include <trianglesystem/gmtrianglesystem>


void mesh::random(float rad, int triangle)
{
    double s=std::sqrt((4*M_PI*rad*rad)/(std::sqrt(3)*triangle));
    //std::cout<<s<<std::endl;
    int k = (M_2PI*rad)/s;
    //std::cout<<k<<std::endl;
    //int edg=3*point-3-k;
    int innerpoint = (triangle+2+k)/2;
    auto vec = GMlib::Vector<float,2>(rad,0);
    GMlib::Angle a=M_2PI/k;
    GMlib::SqMatrix<float,2> matrix (a, GMlib::Vector<float,2>(1.0f,0.0f),GMlib::Vector<float,2>(0.0f,1.0f));

    for (int i=0;i<k;i++)
    {   vec = matrix*vec;
        //std::cout<<vec<<std::endl;
        this->insertAlways(GMlib::TSVertex<float>(static_cast<GMlib::Point<float,2>>(vec)));
    }


    double err = rad-rad*cos(GMlib::Angle(a/2.0).getRad());
    double new_rad = rad-2*err;

    for (int j=0;j<innerpoint;j++)
    {
        auto point = GMlib::TSVertex<float>(gen_float(-new_rad,new_rad),gen_float(-new_rad,new_rad));
        if (checkP(&point, s, new_rad))
            this->insertAlways(point);
        else
            j--;
        //std::cout<<point<<endl;
    }

}

bool mesh::checkP(GMlib::TSVertex<float> *p, float s, float radius)
{
    if (p->getParameter().getLength()>radius)
        return false;

    bool inside = true;

    for (int i = 0; i < this->size(); i++)
    {
        auto op = (*this)[i];

        if ((op.getPos() - p->getPos()).getLength() < (s / 2))
            inside = false;
    }

    return inside;
}


float mesh::gen_float(float a, float b)
{
    float random = (float)std::rand() / (float) RAND_MAX;   //RAND_MAX is a constant defined in <cstdlib>
    float diff = b-a;
    float r = random*diff;
    return a+r;

}

void mesh::regular(int circle, double firstRad, int innerNodes)
{
    this->insertAlways(GMlib::TSVertex<float>(GMlib::Point<float,2>(0,0)));  // 插入中心点

    for (int j=1;j<=circle;j++)
    {
        auto vec = GMlib::Vector<float,2>(firstRad*float(j),0);  // 定义最内圈的半径向量（半径参数是圈数乘以第一圈半径的大小，0度）

        GMlib::Angle a(M_2PI/(innerNodes*j));   //定义一个角度li
        GMlib::SqMatrix<float,2> matrix (a,GMlib::Vector<float,2>(1,0),GMlib::Vector<float,2>(0,1));  //定义一个矩阵，三个参数分别是角度，x轴向量，y轴向量

        for (int i=1;i<=innerNodes*j;i++)
        {
            vec = matrix*vec;    //角度依次旋转
            this->insertAlways(GMlib::TSVertex<float>(static_cast<GMlib::Point<float,2>>(vec)));
            //std::cout<<this->getVertex(i)<<std::endl;
        }
    }
}



void mesh::findElement()
{
    for (int i=0; i<this->size(); i++) //指的mesh里面的所有成员的大小
    {
        if (!(*this)[i].boundary())    //不是所有边界上的成员都放在一个容器里
            //                nodes+=GMlib::TSVertex<float>((*this)[i]);
            nodes+=&(*this)[i];

    }

    A.setDim(nodes.size(),nodes.size());  // 设置矩阵A的大小

    b.setDim(nodes.size());

    for(int i = 0; i < nodes.size(); ++i)
        b[i] = 0.0f;

    for(int i = 0; i < nodes.size(); ++i)
        for(int j = 0; j < nodes.size(); ++j)
            A[i][j] = 0;

    for(int i=0;i<nodes.size();i++)
    {
        for (int j=i+1; j<nodes.size();j++)
        {
            GMlib::TSEdge<float>* commonEdge = nullptr;
            auto edges1 = nodes[i]->getEdges();
            auto edges2 = nodes[j]->getEdges();

            for (int k=0; k< edges1.size();k++)
            {
                for (int m=0; m< edges2.size();m++)
                {
                    if (edges1[k] == edges2[m])
                    {
                        commonEdge = edges1[k];

                    }
                }
            }

            if (commonEdge != nullptr)
            {
                auto getVec=findVector(commonEdge);

                double dd=1/(getVec[0]*getVec[0]);
                double area1=std::abs(getVec[0]^getVec[1]);
                double dh1=dd*(getVec[1]*getVec[0]);
                double h1=dd*area1*area1;

                double area2=std::abs(getVec[2]^getVec[0]);
                double dh2=dd*(getVec[2]*getVec[0]);
                double h2=dd*area2*area2;

                A[i][j]=A[j][i]= ((dh1*(1-dh1)/h1)-dd)*(area1)/2 + ((dh2*(1-dh2)/h2)-dd)*(area2)/2;
            }
        }
    }



    for (int i=0; i<nodes.size();i++)
    {
        auto triangles = nodes[i]->getTriangles();
        for (int k=0; k<triangles.size();k++)
        {
            auto DiaVector = findVectord(nodes[i], triangles[k]);  //对角向量
            A[i][i]+=(DiaVector[2]*DiaVector[2])/(2*std::abs(DiaVector[0]^DiaVector[1]));
        }

    }

    for (int i=0; i<nodes.size();i++)
    {
        auto triangles = nodes[i]->getTriangles();
        for (int k=0; k<triangles.size();k++)
        {
            auto bVector = findVectord(nodes[i], triangles[k]);
          //  std::cout << "bvector: " << bVector << std::endl;
            b[i]+= (std::abs(bVector[0]^bVector[1]))/6;
        }

    }

    //print matrix on screen
    std::cout<<"New matrix"<<std::endl;
    for (int i =0; i<nodes.size();i++)
    {
        for(int j=0; j<nodes.size();j++)
        {
            std::cout<<A[i][j]<<"  ";
        }
        std::cout<<std::endl;
    }

    //print b vector on screen
    std::cout<<"b"<<std::endl;
    for (int i=0; i<nodes.size();i++)
    {

        std::cout<<b[i]<<endl;

    }

    Ainverted=A.invert();
    //std::cout<<Ainverted<<endl;

}

GMlib::Vector<GMlib::Vector<float,2>,3> mesh::findVector(GMlib::TSEdge<float>* commonEdge)   //找三个向量d，a1，a2
{
    GMlib::Vector<GMlib::Vector<float,2>,3> returnVec;


    auto p0=commonEdge->getFirstVertex();      //定义p0是共同边的第一个顶点
    auto p1=commonEdge->getLastVertex();       //定义p1是共同边的最后一个顶点
    GMlib::TSVertex<float> *p2, *p3;           //点定义两个点，是GMlib::TSVertex<float>类型的
    auto triangles=commonEdge->getTriangle();  //定义三角形，是共同边得到的三角形，有两个（库里面的程序）

    auto pointstr1=triangles[0]->getVertices();   //定义一个点，是第一个三角形的顶点
    auto pointstr2=triangles[1]->getVertices();   //定义一个点，是第二个三角的顶点

    for(int i=0;i<3;i++)                          //比较三角形中的三个点
    {
        if(pointstr1[i]!=p0 && pointstr1[i]!=p1)
        {
            p2=pointstr1[i];
        }
        if(pointstr2[i]!=p0 && pointstr2[i]!=p1)
        {
            p3=pointstr2[i];
        }

    }

    GMlib::Vector<float,2> d;
    GMlib::Vector<float,2> a1;
    GMlib::Vector<float,2> a2;

    d=(p1->getParameter())-(p0->getParameter());
    a1=(p2->getParameter())-(p0->getParameter());
    a2=(p3->getParameter())-(p0->getParameter());

    returnVec[0] = d;
    returnVec[1] = a1;
    returnVec[2] = a2;

    return returnVec;

}

GMlib::Vector<GMlib::Vector<float,2>,3> mesh::findVectord(GMlib::TSVertex<float> *point,GMlib::TSTriangle<float>* triangle )
{
    GMlib::Vector<GMlib::Vector<float,2>,3> returnDia;
    auto trianglepoints = triangle->getVertices();

    for(int i=0; i<3;i++)
    {
        if (point == trianglepoints[1])
        {
            std::swap(trianglepoints[0],trianglepoints[1]);
            std::swap(trianglepoints[1],trianglepoints[2]);
        }

        if (point == trianglepoints[2])
        {
            std::swap(trianglepoints[0],trianglepoints[2]);
            std::swap(trianglepoints[1],trianglepoints[2]);
        }


        GMlib::Vector<float,2> d1;
        GMlib::Vector<float,2> d2;
        GMlib::Vector<float,2> d3;


        d1=trianglepoints[2]->getParameter()-trianglepoints[0]->getParameter();
        d2=trianglepoints[1]->getParameter()-trianglepoints[0]->getParameter();
        d3=trianglepoints[2]->getParameter()-trianglepoints[1]->getParameter();

        returnDia[0] = d1;
        returnDia[1] = d2;
        returnDia[2] = d3;

        return returnDia;
    }
}


void mesh::localSimulate(double dt)
{

    if(swp)
    {
        force += dt;
        if(force>1)
            swp=false;
    }
    else
    {
        force -= dt;
        if(force<-1)
            swp=true;
    }

    x= Ainverted*b*force;

    for(int i=0;i<nodes.size();i++)
    {
        nodes[i]->setZ(x[i]);
    }
    this->replot();

}
