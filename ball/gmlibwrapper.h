
//class TestTorus;
class GLContextSurfaceWrapper;
class Ball;
class PhysicalController;
// gmlib
#include <core/gmpoint>

namespace GMlib {

  class Scene;
  class Camera;
  class PointLight;
  class DefaultRenderer;
  class DefaultSelectRenderer;
  class TextureRenderTarget;

  template<typename T, int n>
  class PSurfTexVisualizer;

  template<typename T, int n>
  class PSurf;
}

// qt
#include <QObject>
#include <QSize>
#include <QRectF>
class QOpenGLContext;
class QOpenGLFramebufferObject;
class QMouseEvent;
class QKeyEvent;
class QWheelEvent;


// stl
#include <memory>
#include <unordered_map>


struct RenderCamPair {
  RenderCamPair() : render{nullptr}, camera{nullptr} {}
  std::shared_ptr<GMlib::DefaultRenderer>     render;
  std::shared_ptr<GMlib::Camera>              camera;
  struct {
    QRectF                      geometry { QRectF(0,0,200,200) };
    bool                        changed {true};
  } viewport;
};





class GMlibWrapper : public QObject {
  Q_OBJECT
public:
  explicit GMlibWrapper(std::shared_ptr<GLContextSurfaceWrapper> wrapper );
  ~GMlibWrapper();

  void                                  start();
  void                                  stop();

  const std::shared_ptr<GMlib::Scene>&  getScene() const;
  const GMlib::TextureRenderTarget&     getRenderTextureOf( const std::string& name ) const;

  void                                  initScene();


public slots:
  void                                  changeRenderGeometry( const QString& name,
                                                              const QRectF &new_geometry );



  void                                  mousePressed( const QString& name, QMouseEvent* event );
  void                                  mouseReleased( const QString& name, QMouseEvent* event );
  void                                  mouseDoubleClicked( const QString& name, QMouseEvent* event );
  void                                  mouseMoved( const QString& name, QMouseEvent* event );
  void                                  keyPressed( const QString& name, QKeyEvent* event );
  void                                  keyReleased( const QString& name, QKeyEvent* event );
  void                                  wheelEventOccurred( const QString& name, QWheelEvent* event );

protected:
  void                                  timerEvent(QTimerEvent *e);


private:




  int                                               _timer_id;

  std::shared_ptr<GLContextSurfaceWrapper>          _glsurface;

  std::shared_ptr<GMlib::Scene>                     _scene;
  std::unordered_map<std::string, RenderCamPair>    _rc_pairs;
  std::shared_ptr<GMlib::DefaultSelectRenderer>     _select_renderer;

  std::vector<GMlib::PSurfTexVisualizer<float,3> *>  _textures;


  int                                               _replot_low_medium_high {1};
  bool                                              _move_object_button_pressed {false};
  bool                                              _select_multiple_objects_pressed {false};
  bool                                              _rotate_object_button_pressed {false};
  GMlib::Point<int,2>                               _current_mouse_pos;
  GMlib::Point<int,2>                               _prev_mouse_pos;
  PhysicalController*                               _controller;

  Ball*                                             _contrBall;

signals:
  void                                              signFrameReady();



  // "Singleton instance"
private:
  static std::unique_ptr<GMlibWrapper>      _instance;
public:
  static const GMlibWrapper&                getInstance();
};
