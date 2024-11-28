#ifndef C3T3_INFO_TYPE_H
#define C3T3_INFO_TYPE_H

#include "SMesh_type.h"
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

// Parallel tag
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

typedef CGAL::Polyhedral_mesh_domain_with_features_3<
          EPICK, SMesh, CGAL::Default, int> Polyhedral_mesh_domain;

// Triangulation
typedef CGAL::Mesh_vertex_base_3<EPICK, Polyhedral_mesh_domain,
            CGAL::Triangulation_vertex_base_with_info_3<int, EPICK,
                    CGAL::Regular_triangulation_vertex_base_3<EPICK>>>                               VbInfo;
typedef CGAL::Mesh_cell_base_3<EPICK, Polyhedral_mesh_domain,
            CGAL::Triangulation_cell_base_with_info_3<int, EPICK,
                    CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<EPICK,
                            CGAL::Regular_triangulation_cell_base_3<EPICK>>>>                        CbInfo;
typedef CGAL::Mesh_triangulation_3<Polyhedral_mesh_domain,EPICK,Concurrency_tag,VbInfo,CbInfo>::type TrInfo;
typedef CGAL::Mesh_complex_3_in_triangulation_3<TrInfo>                                              C3t3Info;

typedef TrInfo::Geom_traits Geom_traits;

#endif // C3T3_INFO_TYPE_H
