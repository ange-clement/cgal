#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include "Messages_interface.h"
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

using namespace CGAL::Three;
class CGAL_Lab_triangulate_facets_plugin :
  public QObject,
  public CGAL_Lab_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "triangulate_facets_plugin.json")

  typedef Scene_surface_mesh_item::Face_graph FaceGraph;
  typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

  struct Selection_updater_visitor
    : public CGAL::Polygon_mesh_processing::Hole_filling::Default_visitor
  {
    typedef typename Scene_polyhedron_selection_item::Selection_set_facet Container;
    Container& faces;

    Selection_updater_visitor(Container& container)
      : faces(container)
    {}
    void before_subface_creations(face_descriptor fd)
    {
      Container::iterator it = faces.find(fd);
      faces.erase(it);
    }
    void after_subface_created(face_descriptor fd)
    {
      faces.insert(fd);
    }
    void after_subface_creations() {}
  };

public:

  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;
    actionTriangulateFacets = new QAction("Triangulate Facets", mw);
    actionTriangulateFacets->setProperty("subMenuName","Polygon Mesh Processing");
    if(actionTriangulateFacets) {
      connect(actionTriangulateFacets, SIGNAL(triggered()),
              this, SLOT(triangulate()));
    }
  };

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTriangulateFacets;
  }

  bool applicable(QAction*) const {
    for(CGAL::Three::Scene_interface::Item_id index : scene->selectionIndices()){
      if ( qobject_cast<Scene_surface_mesh_item*>(scene->item(index)) )
        return true;
      if ( qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
        return true;
      if(qobject_cast<Scene_polygon_soup_item*>(scene->item(index)))
        return true;
    }
    return false;
  }

public Q_SLOTS:
   void triangulate()
   {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      for(CGAL::Three::Scene_interface::Item_id index : scene->selectionIndices())
    {
      Scene_surface_mesh_item* sm_item =
          qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

      Scene_polyhedron_selection_item* selection_item =
        qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

      Scene_polygon_soup_item* soup_item =
        qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

      if (soup_item)
      {
        soup_item->triangulate();
      }
      else
      {
        SMesh* pMesh = (sm_item != nullptr) ? sm_item->polyhedron()
                                            : selection_item->polyhedron();

        if(!pMesh)
          continue;

        if(is_triangle_mesh(*pMesh))
        {
          CGAL::Three::Three::warning(tr("The polyhedron  \"%1\"  is already triangulated.")
                            .arg(sm_item->name()) );
          continue;
        }

        if (sm_item)
        {
          if (!CGAL::Polygon_mesh_processing::triangulate_faces(*pMesh))
            CGAL::Three::Three::warning(tr("Some facets could not be triangulated."));

          sm_item->invalidateOpenGLBuffers();
        }
        else if (selection_item)
        {
          Selection_updater_visitor visitor(selection_item->selected_facets);
          if (!CGAL::Polygon_mesh_processing::triangulate_faces(
                  selection_item->selected_facets,
                  *pMesh,
                  CGAL::parameters::visitor(visitor)))
            CGAL::Three::Three::warning(tr("Some facets could not be triangulated."));

          sm_item = selection_item->polyhedron_item();
          selection_item->set_num_faces(num_faces(*sm_item->face_graph()));

          selection_item->invalidateOpenGLBuffers();
          selection_item->itemChanged();
        }

        sm_item->resetColors(); // @todo should have a visitor to give the color of the parent face
      }

      scene->itemChanged(sm_item);
    } // end of the loop on the selected items

    // default cursor
    QApplication::restoreOverrideCursor();
  }

private:
  QAction* actionTriangulateFacets;
  Messages_interface* messages;
};

#include "Triangulate_facets_plugin.moc"
