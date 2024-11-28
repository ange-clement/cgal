#ifndef SCENE_C3T3_COLOR_ITEM_H
#define SCENE_C3T3_COLOR_ITEM_H

#include "Scene_c3t3_color_item_config.h"
#include "C3t3_type.h"

#include "Scene_c3t3_item.h"
// #include "Scene_triangulation_3_item.h"

struct Color_infos
{
    std::vector<float> tetrahedra_colors;
    std::vector<float> facet_colors;
    std::vector<float> edges_colors;
    std::vector<float> vertices_colors;
};

using namespace CGAL::Three;
  class SCENE_C3T3_COLOR_ITEM_EXPORT Scene_c3t3_color_item
  : public Scene_c3t3_item
  {
    Q_OBJECT
  public:
    typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

    Scene_c3t3_color_item(bool is_surface = false) : Scene_c3t3_item(is_surface) { }
    Scene_c3t3_color_item(const C3t3& c3t3, bool is_surface = false) : Scene_c3t3_item(c3t3, is_surface) { }
    virtual ~Scene_c3t3_color_item() { }

    const Color_infos& color_infos() const { return _color_infos; }
    Color_infos& color_infos() { return _color_infos; }

  // public Q_SLOTS:

  protected:

    Color_infos _color_infos;
  };

#endif // SCENE_C3T3_COLOR_ITEM_H
