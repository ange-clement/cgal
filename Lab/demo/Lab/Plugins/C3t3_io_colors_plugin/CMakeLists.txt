include(CGALlab_macros)

# add_library(scene_c3t3_color_item SHARED Scene_c3t3_color_item.cpp)
# target_link_libraries(
#   scene_c3t3_color_item PUBLIC demo_framework CGAL::CGAL_Qt6 Qt6::OpenGLWidgets Qt6::Gui
#                       Qt6::Widgets)
# if(TARGET CGAL::TBB_support)
#     target_link_libraries(scene_c3t3_color_item PUBLIC CGAL::TBB_support)
# endif()

cgal_lab_plugin(c3t3_io_colors_plugin C3t3_io_colors_plugin KEYWORDS
                       Viewer C3t3_io_colors)
# target_link_libraries(c3t3_io_colors_plugin PUBLIC scene_c3t3_color_item)
target_link_libraries(c3t3_io_colors_plugin PUBLIC scene_c3t3_item)

# if(TBB_FOUND)
#   target_link_libraries(c3t3_io_colors_plugin PUBLIC CGAL::TBB_support)
# endif()
# target_link_libraries(
#   phd_plugin PUBLIC scene_polygon_soup_item scene_image_item )
