#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = typename Kernel::Point_3;
using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;



template <typename Domain>
class Label_cell_offset
{
  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Vertex_descriptor = typename Domain::Vertex_descriptor;

public:
  Label_cell_offset(Vector_3 offset)
      : offset(offset)
  { }

  bool position(const Domain& domain,
                const typename Domain::Cell_descriptor& vh,
                typename Domain::Geom_traits::Point_3& point)
  {
    typename Domain::Cell_vertices vertices = domain.cell_vertices(vh);

    point = domain.point(vertices[0]) + offset;

    FT v0 = domain.value(vertices[0]);
    for (unsigned int i = 1, size = vertices.size(); i < size; i++)
    {
      if (domain.value(vertices[i]) != v0) {
        return true;
      }
    }
    return false;
  }

private:
  Vector_3 offset;
};

class Label_condition
{
public:
  template <typename Domain,
            typename Edge_descriptor>
  bool condition(const Domain& domain, const Edge_descriptor & e, bool & orientation)
  {
    using FT = typename Domain::Geom_traits::FT;
    const auto& vertices = domain.incident_vertices(e);
    const FT s0 = domain.value(vertices[0]);
    const FT s1 = domain.value(vertices[1]);

    if (s0 == s1)
    {
      return false;
    }

    orientation = s0 < s1;
    return true;
  }
};


int main(int, char**)
{
  const std::string fname = CGAL::data_file_path("images/liver.inr.gz");

  // load volumetric image from a file
  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read image file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // convert image to a Cartesian grid
  Grid grid{image};

  // create a domain from the grid
  using Domain = CGAL::Isosurfacing::Explicit_Cartesian_grid_domain_3<Grid>;
  Domain domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  // prepare collections for the output indexed mesh
  Point_range points;
  Polygon_range polygons;

  // create condition and position
  Kernel::Vector_3 offset(.5*image.vx(), .5*image.vy(), .5*image.vz());
  Label_cell_offset<Domain> position(offset);
  Label_condition condition;

  // execute dual_contouring
  CGAL::Isosurfacing::dual_contouring(domain, points, polygons, condition, position);

  // save output indexed mesh to a file, in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
