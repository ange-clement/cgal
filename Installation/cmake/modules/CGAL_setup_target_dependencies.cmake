function(CGAL_setup_target_dependencies target)
  if(${target} STREQUAL CGAL)
    CGAL_setup_CGAL_dependencies(${target})
  elseif(${target} STREQUAL CGAL_Core)
    CGAL_setup_CGAL_Core_dependencies(${target})
  elseif(${target} STREQUAL CGAL_ImageIO)
    CGAL_setup_CGAL_ImageIO_dependencies(${target})
  elseif(${target} STREQUAL CGAL_Qt6)
    CGAL_setup_CGAL_Qt6_dependencies(${target})
  endif()
endfunction()
