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


class Label_cell_center
{
public:
  template <typename Domain>
  bool position(const Domain& domain,
                const typename Domain::Cell_descriptor& vh,
                typename Domain::Geom_traits::Point_3& point) const
  {
    using Geom_traits = typename Domain::Geom_traits;
    using FT = typename Geom_traits::FT;
    using Point_3 = typename Geom_traits::Point_3;
    using Vector_3 = typename Geom_traits::Vector_3;

    using Vertex_descriptor = typename Domain::Vertex_descriptor;

    typename Domain::Cell_vertices vertices = domain.cell_vertices(vh);

    std::vector<Point_3> pos(vertices.size());
    std::transform(vertices.begin(), vertices.end(), pos.begin(),
                   [&](const Vertex_descriptor& v) { return domain.point(v); });

    // set point to cell center
    point = CGAL::centroid(pos.begin(), pos.end(), CGAL::Dimension_tag<0>());

    FT v0 = domain.value(vertices[0]);
    for (unsigned int i = 1, size = vertices.size(); i < size; i++)
    {
      if (domain.value(vertices[i]) != v0) {
        return true;
      }
    }
    return false;
  }
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
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain(grid);

  // prepare collections for the output indexed mesh
  Point_range points;
  Polygon_range polygons;

  // create condition and position
  Label_cell_center position;
  Label_condition condition;

  // execute dual_contouring
  CGAL::Isosurfacing::dual_contouring(domain, points, polygons, condition, position);

  // save output indexed mesh to a file, in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
