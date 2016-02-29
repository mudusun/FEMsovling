#include "gmlibwrapper.h"
#include "glcontextsurfacewrapper.h"

#include "plane.h"
#include "ball.h"
#include "physicalcontroller.h"
#include "collisionobject.h"
//#include "testtorus.h"
#include "utils.h"


// GMlib
#include <gmOpenglModule>
#include <gmSceneModule>
#include <gmParametricsModule>


// Qt
#include <QTimerEvent>
#include <QRectF>
#include <QMouseEvent>
#include <QDebug>
#include <QImage>

#include <QImageReader>
#include <QDirIterator>

// stl
#include <stdexcept>
#include <thread>
#include <mutex>





std::unique_ptr<GMlibWrapper> GMlibWrapper::_instance {nullptr};


GMlibWrapper::GMlibWrapper(std::shared_ptr<GLContextSurfaceWrapper> context)
//  : GMlibWrapper()
  : QObject(), _timer_id{0}, _glsurface(context), _select_renderer{nullptr}
{

  if(_instance != nullptr) {

    std::cerr << "This version of the GMlibWrapper only supports a single instance of the GMlibWraper..." << std::endl;
    std::cerr << "Only one of me(0x" << this << ") please!! ^^" << std::endl;
    assert(!_instance);
    exit(666);
  }


  _instance = std::unique_ptr<GMlibWrapper>(this);

  _glsurface->makeCurrent();

  // Setup and initialized GMlib GL backend
  GMlib::GL::OpenGLManager::init();

  // Setup and init the GMlib GMWindow
  _scene = std::make_shared<GMlib::Scene>();
}

GMlibWrapper::~GMlibWrapper() {

  stop();

  _glsurface->makeCurrent(); {

    _select_renderer->releaseCamera();
    _select_renderer.reset();

    for( auto& rc_pair : _rc_pairs ) {

      rc_pair.second.render->releaseCamera();
      _scene->removeCamera( rc_pair.second.camera.get() );

      rc_pair.second.render.reset();
      rc_pair.second.camera.reset();
    }

    _scene->clear();

  } _glsurface->doneCurrent();
}

void GMlibWrapper::changeRenderGeometry(const QString& name, const QRectF& geometry) {

  const QSize size = geometry.size().toSize();

  if( size.width() <= 0 || size.height() <= 0 )
    return;

  if( _rc_pairs.count(name.toStdString()) <= 0 )
    return;

  auto& rc_pair = _rc_pairs[name.toStdString()];
  if(rc_pair.viewport.geometry == geometry )
    return;

  rc_pair.viewport.geometry = geometry;
  rc_pair.viewport.changed = true;
}

void GMlibWrapper::timerEvent(QTimerEvent* e) {

  e->accept();

  // Simuation order
  // 1) Prepare must be run first
  // 2) Simulate and render can be run in parallell


  // Grab and activate GL context
  _glsurface->makeCurrent(); {

    // 1)
    _scene->prepare();


    _scene->simulate();

//    std::vector<std::thread> threads;

    // Add simulation thread
//    threads.push_back(std::thread(&GMlib::Scene::simulate,_scene));

    // Add Render threads
    for( auto& rc_pair : _rc_pairs ) {
  //      qDebug() << "About to render: " << rc_pair.first.c_str();
  //      qDebug() << "  Viewport: ";
  //      qDebug() << "    Changed: " << rc_pair.second.viewport.changed;
  //      qDebug() << "    Geometry: " << rc_pair.second.viewport.geometry;

      if(rc_pair.second.viewport.changed) {
        const QSizeF size = rc_pair.second.viewport.geometry.size();
        rc_pair.second.render->reshape( GMlib::Vector<int,2>(size.width(),size.height()));
        rc_pair.second.camera->reshape( 0, 0, size.width(), size.height() );
        rc_pair.second.viewport.changed = false;
      }

      rc_pair.second.render->render();
      rc_pair.second.render->swap();
    }

//    for( auto& thread : threads )
//      thread.join();

  } _glsurface->doneCurrent();

  emit signFrameReady();
}

const GMlibWrapper&
GMlibWrapper::getInstance() {
  return *_instance;
}

void GMlibWrapper::start() {

  if( _timer_id || _scene->isRunning() )
    return;

  _timer_id = startTimer(16, Qt::PreciseTimer);
  _scene->start();
}

void GMlibWrapper::stop() {

  if( !_timer_id || !_scene->isRunning() )
    return;

  _scene->stop();
  killTimer(_timer_id);
  _timer_id = 0;
}

void GMlibWrapper::initScene() {

  // Make OpenGL context current on offscreensurface
  _glsurface->makeCurrent(); {

    // Insert a light
    GMlib::Point<GLfloat,3> init_light_pos( 2.0, 4.0, 10 );
    GMlib::PointLight *light = new GMlib::PointLight(  GMlib::GMcolor::White, GMlib::GMcolor::White,
                                                       GMlib::GMcolor::White, init_light_pos );
    light->setAttenuation(0.8, 0.002, 0.0008);
    _scene->insertLight( light, false );

    // Insert Sun
    _scene->insertSun();


    int init_viewport_size = 600;
    GMlib::Point<float,3> init_cam_pos(  0.0f, 0.0f, 0.0f );
    GMlib::Vector<float,3> init_cam_dir( 0.0f, 1.0f, 0.0f );
    GMlib::Vector<float,3> init_cam_up(  0.0f, 0.0f, 1.0f );

    _rc_pairs.reserve(4);
    _rc_pairs["Projection"] = RenderCamPair {};
    _rc_pairs["Front"]      = RenderCamPair {};
    _rc_pairs["Side"]       = RenderCamPair {};
    _rc_pairs["Top"]        = RenderCamPair {};

    for( auto& rcpair : _rc_pairs ) {

      rcpair.second.render = std::make_shared<GMlib::DefaultRenderer>();
      rcpair.second.camera = std::make_shared<GMlib::Camera>();
      rcpair.second.render->setCamera(rcpair.second.camera.get());
    }

    // Projection cam
    auto& proj_rcpair = _rc_pairs["Projection"];
    proj_rcpair.camera->set(init_cam_pos,init_cam_dir,init_cam_up);
    proj_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
    proj_rcpair.camera->rotateGlobal( GMlib::Angle(-45), GMlib::Vector<float,3>( 1.0f, 0.0f, 0.0f ) );
    proj_rcpair.camera->translateGlobal( GMlib::Vector<float,3>( 0.0f, -20.0f, 20.0f ) );
    _scene->insertCamera( proj_rcpair.camera.get() );
    proj_rcpair.render->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );
    proj_rcpair.render->setClearColor(GMlib::GMcolor::Black);

    // Front cam
    auto& front_rcpair = _rc_pairs["Front"];
    front_rcpair.camera->set( init_cam_pos + GMlib::Vector<float,3>( 0.0f, -50.0f, 0.0f ), init_cam_dir, init_cam_up );
    front_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
    _scene->insertCamera( front_rcpair.camera.get() );
    front_rcpair.render->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );

    // Side cam
    auto& side_rcpair = _rc_pairs["Side"];
    side_rcpair.camera->set( init_cam_pos + GMlib::Vector<float,3>( -50.0f, 0.0f, 0.0f ), GMlib::Vector<float,3>( 1.0f, 0.0f, 0.0f ), init_cam_up );
    side_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
    _scene->insertCamera( side_rcpair.camera.get() );
    side_rcpair.render->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );

    // Top cam
    auto& top_rcpair = _rc_pairs["Top"];
    top_rcpair.camera->set( init_cam_pos + GMlib::Vector<float,3>( 0.0f, 0.0f, 50.0f ), -init_cam_up, init_cam_dir );
    top_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
    _scene->insertCamera( top_rcpair.camera.get() );
    top_rcpair.render->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );





//    // Iso Camera
//    auto& isorcpair = (_rc_pairs["Iso"] = RenderCamPair {});
//    isorcpair.render = std::make_shared<GMlib::DefaultRenderer>();
//    isorcpair.camera = std::make_shared<GMlib::IsoCamera>();
//    isorcpair.render->setCamera(isorcpair.camera.get());
//    _scene->insertCamera( isorcpair.camera.get() );
//    isorcpair.camera->set(init_cam_pos,init_cam_dir,init_cam_up);
//    isorcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
//    isorcpair.camera->rotateGlobal( GMlib::Angle(-45), GMlib::Vector<float,3>( 1.0f, 0.0f, 0.0f ) );
//    isorcpair.camera->translate( GMlib::Vector<float,3>( 0.0f, -20.0f, 20.0f ) );
//    isorcpair.render->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );


    // Setup Select Renderer
    _select_renderer = std::make_shared<GMlib::DefaultSelectRenderer>();



//#define TEST_CURVE
#define TEST_SURFACE


#ifdef TEST_CURVE

//    // Curve visualizers
//    auto curve_visualizer = new GMlib::PCurveDerivativesVisualizer<float,3>;

//    // Curve
//    auto curve = new GMlib::PCircle<float>(2.0f);
//    curve->toggleDefaultVisualizer();
//    curve->insertVisualizer(curve_visualizer);
//    curve->replot(100,1);
//    _scene->insert(curve);

#endif





#ifdef TEST_SURFACE

    // Surface visualizers
//    auto surface_visualizer = new GMlib::PSurfDerivativesVisualizer<float,3>;
//    auto surface_visualizer = new GMlib::PSurfNormalsVisualizer<float,3>;
    auto surface_visualizer = new GMlib::PSurfParamLinesVisualizer<float,3>;
   // auto surface_visualizer = new GMlib::PSurfPointsVisualizer<float,3>;

    // Surface
//    auto surface = new TestTorus;
//    surface->toggleDefaultVisualizer();
//    surface->insertVisualizer(surface_visualizer);
//    surface->replot(200,200,1,1);
//    _scene->insert(surface);

//    surface->test01();


    // create the wall here
    GMlib::Point<float,3> point1(-6,-6,0);
    GMlib::Point<float,3> point2(6,-6,0);
    GMlib::Point<float,3> point3(-6,6,0);

    GMlib::Vector<float,3> vec1(12,0,0);
    GMlib::Vector<float,3> vec2(0,12,0);
    GMlib::Vector<float,3> vec3(0,0,4);   //墙的高度是3


    GMlib::DMatrix<GMlib::Vector<float,3> > m0(4,4);
    GMlib::Vector<float,3> vec000(-6,-6,0); //曲面左下角点的坐标
    GMlib::Vector<float,3> vec001(-2,-6,0);
    GMlib::Vector<float,3> vec002(2,-6,0);
    GMlib::Vector<float,3> vec003(6,-6,0);


    GMlib::Vector<float,3> vec010(-6,-2,0);
    GMlib::Vector<float,3> vec011(-2,-2,-5); // the z in this vector said that the surface convex or concave z轴的正负号表示凸面还是曲面
    GMlib::Vector<float,3> vec012(2,-2,3);
    GMlib::Vector<float,3> vec013(6,-2,0);


    GMlib::Vector<float,3> vec020(-6,2,0);
    GMlib::Vector<float,3> vec021(-2,2,-5);
    GMlib::Vector<float,3> vec022(2,2,3); //右边
    GMlib::Vector<float,3> vec023(6,2,0);

    GMlib::Vector<float,3> vec030(-6,6,0);
    GMlib::Vector<float,3> vec031(-2,6,0);
    GMlib::Vector<float,3> vec032(2,6,0);
    GMlib::Vector<float,3> vec033(6,6,0); //曲面右上角点的坐标

    m0[0][0] = vec000; m0[0][1] = vec001; m0[0][2] = vec002; m0[0][3] = vec003;
    m0[1][0] = vec010; m0[1][1] = vec011; m0[1][2] = vec012; m0[1][3] = vec013;
    m0[2][0] = vec020; m0[2][1] = vec021; m0[2][2] = vec022; m0[2][3] = vec023;
    m0[3][0] = vec030; m0[3][1] = vec031; m0[3][2] = vec032; m0[3][3] = vec033;



//************************  change visualization of wall    ************************************************
   using namespace GMlib;


  //   PSurfNormalsVisualizer<float,3> *normals_visu = new PSurfNormalsVisualizer<float,3>;

  //  QImage img = QImage("<abs_path_to_texture>/chess_oak.jpg");

    for(int i=0; i<4; i++)
    {
    QString texHandle1 = "C:/Users/hxu001/Documents/ball/Image/Brick_wall.jpg";
    QImage img1 = QImage(texHandle1);

    GMlib::PSurfTexVisualizer<float,3> *tex_visu1= new GMlib::PSurfTexVisualizer<float,3>;


    if(! QFile::exists(texHandle1))

        std::cerr << "File " << texHandle1.toStdString() <<" does not exists"<<std::endl;

    if( img1.isNull() )
    {

      qDebug() << QImageReader::supportedImageFormats();

    }

    assert(!img1.isNull());
   // img1 = img1.convertToFormat(QImage::Format_RGB888);


   // const QImage gl_img = QGLWidget::convertToGLFormat( img1 );


    GMlib::GL::Texture tex1;

    tex1.create( GL_TEXTURE_2D );

    tex1.texImage2D( 0, GL_RGB32F, img1.width(), img1.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, img1.bits() );

    tex1.texParameteri( GL_TEXTURE_MIN_FILTER, GL_LINEAR );

    tex1.texParameteri( GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    tex1.texParameterf( GL_TEXTURE_WRAP_S, GL_REPEAT);

    tex1.texParameterf( GL_TEXTURE_WRAP_T, GL_REPEAT);

    tex_visu1->setTexture(tex1);

    _textures.push_back(tex_visu1);
}


//**************************************************************************************************************

    auto floor = new Plane<float> (m0);
//    auto floor = new GMlib::PBezierSurf<float> (m0);

    floor->toggleDefaultVisualizer();
    floor->replot(20,20,1,1);
    floor->getDefaultVisualizer()->setDisplayMode(GMlib::Visualizer::DISPLAY_MODE_WIREFRAME);
    floor->enableDefaultVisualizer();
    auto ctrl=new PhysicalController(floor);
    _scene->insert(ctrl);


    //plane North
    GMlib::PPlane<float> *plane1 = new GMlib::PPlane<float>(point3,vec1,vec3);
    plane1->insertVisualizer(_textures[0]);
    //plane1->toggleDefaultVisualizer();
    plane1->replot(20,20,1,1);
    plane1->setMaterial(GMlib::GMmaterial::Silver);
    ctrl->addplane(plane1);


    //plane West
    GMlib::PPlane<float> *plane2 = new GMlib::PPlane<float>(point1,vec2,vec3);
    plane2->insertVisualizer(_textures[1]);
    // plane2->toggleDefaultVisualizer();
    plane2->replot(200,200,1,1);
    plane2->setMaterial(GMlib::GMmaterial::Silver);
    ctrl->addplane(plane2);


    //plane East
    GMlib::PPlane<float> *plane3 = new GMlib::PPlane<float>(point2,vec3,vec2);
    plane3->insertVisualizer(_textures[2]);
    //plane3->toggleDefaultVisualizer();
    plane3->replot(20,20,1,1);
    plane3->setMaterial(GMlib::GMmaterial::Silver);
    ctrl->addplane(plane3);



    //plane South
    GMlib::PPlane<float> *plane4 = new GMlib::PPlane<float>(point1,vec3,vec1);
    plane4->insertVisualizer(_textures[3]);
    //plane4->toggleDefaultVisualizer();
    plane4->replot(20,20,1,1);
    plane4->setMaterial(GMlib::GMmaterial::Silver);
    ctrl->addplane(plane4);


    //Add balls
    auto ball1=new Ball(GMlib::Vector<float,3>(-2.0f,6.0f,0.0f),1,1,floor);        //定义球的速度，质量，半径，曲面
    ctrl->addball(ball1);
    ball1->enableDefaultVisualizer();
    ball1->translate(GMlib::Vector<float,3> (-2,-2,0));
    ball1->setMaterial(GMlib::GMmaterial::Emerald);
    ball1->replot(20,20,1,1);


    auto ball2=new Ball(GMlib::Vector<float,3>(2.0f,12.0f,0.0f),1,0.7,floor);        //定义球的速度，质量，半径，曲面
    ctrl->addball(ball2);
    ball2->enableDefaultVisualizer();
    ball2->translate(GMlib::Vector<float,3> (-3,0,0));
    ball2->setMaterial(GMlib::GMmaterial::PolishedRed);
    ball2->replot(20,20,1,1);

    auto ball3=new Ball(GMlib::Vector<float,3>(0.3f,4.0f,0.0f),1,1.2,floor);        //定义球的速度，质量，半径，曲面
    ctrl->addball(ball3);
    ball3->enableDefaultVisualizer();
    ball3->translate(GMlib::Vector<float,3> (0,-2,0));
    ball3->setMaterial(GMlib::GMmaterial::Turquoise);
    ball3->replot(20,20,1,1);

    auto ball4=new Ball(GMlib::Vector<float,3>(0.0f,4.0f,0.0f),1,0.8,floor);        //定义球的速度，质量，半径，曲面
    ctrl->addball(ball4);
    ball4->enableDefaultVisualizer();
    ball4->translate(GMlib::Vector<float,3> (-4,1,0));
    ball4->setMaterial(GMlib::GMmaterial::Pearl);
    ball4->replot(20,20,1,1);

    _contrBall = new Ball(GMlib::Vector<float,3>(0.0f,4.0f,0.0f),1,0.5,floor);
    ctrl->addball(_contrBall);
    _contrBall->insertVisualizer(surface_visualizer);
    _contrBall->replot(100,100,1,1);
    _contrBall->setMaterial(GMlib::GMmaterial::Gold);
    _contrBall->translate(GMlib::Point<float,3>(0,5,1));

#endif

//********************************************************************************

  } _glsurface->doneCurrent();
}

const std::shared_ptr<GMlib::Scene>&
GMlibWrapper::getScene() const {

  return _scene;
}

const GMlib::TextureRenderTarget&
GMlibWrapper::getRenderTextureOf(const std::string& name) const {

  if(!_rc_pairs.count(name)) throw std::invalid_argument("[][]Render/Camera pair '" + name + "'  does not exist!");

  return _rc_pairs.at(name).render->getFrontRenderTarget();
}




///////////////////////////////////////////////////////////////
/// \brief GMlibWrapper::mousePressed
/// \param name  - the name of the camera
/// \param event - the mouse event, using position and butten
///////////////////////////////////////////////////////////////


void GMlibWrapper::mousePressed(const QString& name, QMouseEvent* event )
{
    const QPointF& pos = event->pos();

    _current_mouse_pos = _prev_mouse_pos = GMlib::Point<int,2>(pos.x(),pos.y());

    const auto& rc_select = _rc_pairs.at(name.toStdString());
    const auto& rc_geo    =  rc_select.viewport.geometry;

    auto& rc   = _rc_pairs.at(name.toStdString());
    auto cam   =  rc.camera.get();

    if(event->button() == Qt::RightButton)
    {
      if(_rotate_object_button_pressed)         // Ctrl
      {
        GMlib::SceneObject * sel_obj = getScene()->getSelectedObjects()(0);

        qDebug() << "Ctrl pressed" << sel_obj;

        if(sel_obj)
        {
          cam->lock(sel_obj);
        }
        else
        {
          if(cam->isLocked()) cam->unLock();

          cam->lock((getScene()->getSphereClean().getPos() - cam->getPos())*cam->getDir());

          qDebug() << "Lock" << (getScene()->getSphereClean().getPos() - cam->getPos())*cam->getDir();
        }
      }
      else if(_move_object_button_pressed)      // Shift
      {
        qDebug() << "Shift pressed";
      }
      else if(_select_multiple_objects_pressed) // Alt
      {
        qDebug() << "Alt pressed";
        GMlib::SceneObject* obj{nullptr};

        GMlib::Vector<int,2> size(rc_geo.width(),rc_geo.height());
        _select_renderer->setCamera(cam);

        _glsurface->makeCurrent();
        {
          _select_renderer->reshape(size);
          _select_renderer->prepare();
          _select_renderer->select(0);
          obj = _select_renderer->findObject(pos.x(),size(1)-pos.y()-1);
        }
        _glsurface->doneCurrent();

        if(obj)
          obj->toggleSelected();
      }
      else
      {
        GMlib::SceneObject* obj{nullptr};
        _select_renderer->setCamera(cam);

        GMlib::Vector<int,2> size(rc_geo.width(),rc_geo.height());

        _glsurface->makeCurrent();
        {
          _select_renderer->reshape(size);
          _select_renderer->prepare();
          _select_renderer->select(0);
          obj = _select_renderer->findObject(pos.x(),size(1)-pos.y()-1);
        }
        _glsurface->doneCurrent();

        if(obj)
        {
         // if(obj == cam->getLockObj())
            cam->unLock();
          //else
            cam->lock(obj);
        }
        else
        {
          if(cam->isLocked())
            cam->unLock();
          else
            cam->lock((getScene()->getSphereClean().getPos() - cam->getPos())*cam->getDir());
        }
      }
    }
    qDebug() << getScene()->getSelectedObjects().getSize();
}

void GMlibWrapper::mouseReleased(const QString& name, QMouseEvent* event) {
    Q_UNUSED(name)
    Q_UNUSED(event)
}

void GMlibWrapper::mouseDoubleClicked(const QString& name, QMouseEvent* event) {
    Q_UNUSED(name)
    Q_UNUSED(event)
}

///////////////////////////////////////////////////////////////
/// \brief GMlibWrapper::mouseMoved
/// \param name  - the name of the camera
/// \param event - the mouse event, using position and butten
///////////////////////////////////////////////////////////////

void GMlibWrapper::mouseMoved(const QString &name, QMouseEvent *event)
{
    Q_UNUSED(name);
    Q_UNUSED(event);

    const QPointF& pos = event->pos();

    _current_mouse_pos = GMlib::Point<int,2>(pos.x(), pos.y());

    if(_prev_mouse_pos.getLength() <= 0.0f)
        _prev_mouse_pos = _current_mouse_pos;

    auto& rc    = _rc_pairs.at(name.toStdString());
    auto  cam   = rc.camera.get();
    float SNAP  = 0.01f;

    if(_move_object_button_pressed)
    {
        _glsurface->makeCurrent();

        const GMlib::Array<GMlib::SceneObject*> &sel_objs = getScene()->getSelectedObjects();

        for( int i = 0; i < sel_objs.getSize(); i++ ) {
            GMlib::SceneObject* obj = sel_objs(i);

            if( obj ) {
                const double dh = cam->deltaTranslate( obj );
                const GMlib::Vector<float,3> deltav(
                            ( ( _prev_mouse_pos(0) - _current_mouse_pos(0) ) * dh ) * cam->getSide() +
                            ( ( _prev_mouse_pos(1) - _current_mouse_pos(1) ) * dh ) * cam->getUp() );

                if( deltav.getLength() > SNAP && deltav.getLength() < 1000.0 ) {
                    if( obj->getTypeId() != GMlib::GM_SO_TYPE_SELECTOR )
                        obj->translateGlobal( deltav );
                    else if( obj->getTypeId()== GMlib::GM_SO_TYPE_SELECTOR )
                        obj->editPos(deltav);
                }
            }
        }
        _glsurface->doneCurrent();
    }
    else if(_rotate_object_button_pressed)
    {
        _glsurface->makeCurrent();

        const GMlib::Vector<int,2> pos = _current_mouse_pos;
        const GMlib::Vector<int,2> prev = _prev_mouse_pos;
        const GMlib::Array<GMlib::SceneObject*> &objs = getScene()->getSelectedObjects();

        // Compute rotation axis and angle in respect to the camera and view.
        const GMlib::UnitVector<float,3> rot_v =
                float( pos(0) - prev(0) ) * cam->getUp() -
                float( pos(1) - prev(1) ) * cam->getSide();

        const GMlib::Angle ang(
                M_2PI * sqrt(
                   pow( double( pos(0) - prev(0) ) / cam->getViewportW(), 2 ) +
                   pow( double( prev(1) - pos(1) ) / cam->getViewportH(), 2 ) ) );

        int no_objs = 0;
        GMlib::Sphere<float,3> sphere;
        for( int i = 0; i < objs.getSize(); ++i )
            if( objs(i)->getTypeId() != GMlib::GM_SO_TYPE_SELECTOR ) {
                sphere += objs(i)->getSurroundingSphereClean();
                no_objs++;
            }

        for( int i = 0; i < objs.getSize(); ++i )

            if( objs(i)->getTypeId() != GMlib::GM_SO_TYPE_SELECTOR )
                if( std::abs(pos(0)-prev(0)) > POS_TOLERANCE || std::abs(pos(1)-prev(1)) > POS_TOLERANCE )
                    no_objs > 1 ? objs(i)->rotate( ang, sphere.getPos(), rot_v) : objs(i)->rotate( ang, rot_v);

        _glsurface->doneCurrent();
    }
    else if(_select_multiple_objects_pressed) // Alt
    {

    }
    else
    {
      float scale;

      if(cam->isLocked()) scale = M_2PI * cam->getLockDist();
      else                scale = getScene()->getSphere().getRadius();

      const GMlib::Vector<float,2> delta (
            (_current_mouse_pos(0) - _prev_mouse_pos(0)) * scale / cam->getViewportW(),
            (_current_mouse_pos(1) - _prev_mouse_pos(1)) * scale / cam->getViewportH() );

      cam->move( delta );
    }

    _prev_mouse_pos = _current_mouse_pos;
}
void
GMlibWrapper::keyPressed(const QString& name, QKeyEvent* event) {
    Q_UNUSED(name)

    if( event->key() == Qt::Key_R ) _scene->toggleRun();

    if( event->key() == Qt::Key_E )
    {
        _glsurface->makeCurrent();

        const GMlib::Array<GMlib::SceneObject*> &sel_objs = getScene()->getSelectedObjects();
        for( int i = 0; i < sel_objs.getSize(); i++ ) {

            GMlib::SceneObject *sel_obj = sel_objs(i);

            // ERBS
            GMlib::PERBSSurf<float> *esObj = dynamic_cast<GMlib::PERBSSurf<float>*>( sel_obj );

            // Bezier
            GMlib::PBezierSurf<float> *bsObj = dynamic_cast<GMlib::PBezierSurf<float>*>( sel_obj );


            if( esObj ) {

                if( esObj->isLocalPatchesVisible() )
                    esObj->hideLocalPatches();
                else
                    esObj->showLocalPatches();

            }
            else if( bsObj ) {

                GMlib::PERBSSurf<float> *parent = dynamic_cast<GMlib::PERBSSurf<float>*>( bsObj->getParent());
                if( parent ) {

                    if( bsObj->toggleCollapsed() )
                        bsObj->hideSelectors();
                    else
                        bsObj->showSelectors(true);
                }
                else {

                    if( bsObj->isSelectorsVisible() )
                        bsObj->hideSelectors();
                    else
                        bsObj->showSelectors(true);
                }

            }
        }
        _glsurface->doneCurrent();
    }

    if(event->key() == Qt::Key_Shift)
    {
        _move_object_button_pressed = true;
    }

    if(event->key() == Qt::Key_Control)
    {
        _rotate_object_button_pressed = true;
    }

    if(event->key() == Qt::Key_Alt)
    {
        _select_multiple_objects_pressed = true;
    }

    if(event->key() == Qt::Key_W)
    {
        _glsurface->makeCurrent();

        const GMlib::Array<GMlib::SceneObject*> &sel_objs = getScene()->getSelectedObjects();
        for( int i = 0; i < sel_objs.getSize(); i++ ) {

            GMlib::SceneObject *sel_obj = sel_objs(i);

            // ERBS
            GMlib::PERBSSurf<float> *es_obj = dynamic_cast<GMlib::PERBSSurf<float>*>( sel_obj );

            if( es_obj ) {

              GMlib::Array<GMlib::Visualizer*> &visus = es_obj->getVisualizers();
              for( int i = 0; i < visus.getSize(); ++i ) {
                visus[i]->toggleDisplayMode();
              }
              es_obj->replot();
            }
        }
        _glsurface->doneCurrent();
    }

    if(event->key() == Qt::Key_P)
    {
        _glsurface->makeCurrent();

        const GMlib::Array<GMlib::SceneObject*> &sel_objs = getScene()->getSelectedObjects();
        for( int i = 0; i < sel_objs.getSize(); i++ ) {

            GMlib::SceneObject *sel_obj = sel_objs(i);

            // ERBS
            GMlib::PERBSSurf<float> *es_obj = dynamic_cast<GMlib::PERBSSurf<float>*>( sel_obj );
            _replot_low_medium_high += 2;

            if(_replot_low_medium_high>5)
                _replot_low_medium_high = 1;

            if( es_obj )
                es_obj->replot(_replot_low_medium_high*10,_replot_low_medium_high*10,2,2);
        }
        _glsurface->doneCurrent();
    }

   //***************************************************
      if(event->key() == Qt::Key_Up)
        {
            _contrBall->moveUp();
        }

        if(event->key() == Qt::Key_Down)
        {
            _contrBall->moveDown();
        }

        if(event->key() == Qt::Key_Left)
        {
            _contrBall->moveLeft();
        }

        if(event->key() == Qt::Key_Right)
        {
            _contrBall->moveRight();
        }
  //****************************************************************************

}

void GMlibWrapper::keyReleased(const QString& name, QKeyEvent* event) {
    Q_UNUSED(name)

    if(event->key() == Qt::Key_Shift)
        _move_object_button_pressed = false;
    if(event->key() == Qt::Key_Control)
        _rotate_object_button_pressed = false;
    if(event->key() == Qt::Key_Alt)
        _select_multiple_objects_pressed = false;
}

void GMlibWrapper::wheelEventOccurred(const QString& name, QWheelEvent* event) {

    int delta = event->delta();
    auto& rc = _rc_pairs.at(name.toStdString());

    const auto& camera_geo = rc.viewport.geometry;
    auto camera    = rc.camera.get();
    auto isocamera = dynamic_cast<GMlib::IsoCamera*>(camera);
    if(isocamera) {

        if( delta < 0 ) isocamera->zoom( 1.05 );
        if( delta > 0 ) isocamera->zoom( 0.95 );
    }
    else if(camera) {

        double scale;
        if( camera->isLocked() )
            scale = camera->getLockDist();
        else
            scale = getScene()->getSphere().getRadius();

        camera->move( delta*scale / camera_geo.height());

    }

}
