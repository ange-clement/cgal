#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <memory>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;
typedef CGAL::Straight_skeleton_2<K>  Ss ;

typedef std::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon_2 outer ;
  outer.push_back( Point(-1,-1) ) ;
  outer.push_back( Point(0,-12) ) ;
  outer.push_back( Point(1,-1) ) ;
  outer.push_back( Point(12,0) ) ;
  outer.push_back( Point(1,1) ) ;
  outer.push_back( Point(0,12) ) ;
  outer.push_back( Point(-1,1) ) ;
  outer.push_back( Point(-12,0) ) ;

  Polygon_2 hole ;
  hole.push_back( Point(-1,0) ) ;
  hole.push_back( Point(0,1 ) ) ;
  hole.push_back( Point(1,0 ) ) ;
  hole.push_back( Point(0,-1) ) ;

  assert(outer.is_counterclockwise_oriented());
  assert(hole.is_clockwise_oriented());

  Polygon_with_holes poly( outer ) ;
  poly.add_hole( hole ) ;

  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly);

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);

  return EXIT_SUCCESS;
}
