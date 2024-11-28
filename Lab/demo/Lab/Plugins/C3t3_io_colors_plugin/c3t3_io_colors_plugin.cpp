#include <CGAL/SMDS_3/io_signature.h>
// #include "Scene_c3t3_color_item.h"
#include "Scene_c3t3_item.h"
#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/IO/File_avizo.h>
#include <CGAL/Three/Triangle_container.h>
#include <iostream>
#include <fstream>

#include <QMessageBox>

struct Color_infos
{
    std::vector<float> vertices_colors;
    std::vector<float> edges_colors;
    std::vector<float> facets_colors;
    std::vector<float> cells_colors;
};

struct Draw_Color_infos
{
    bool flat_rendering_mode;
    std::vector<float> flat_facets_colors;
    std::vector<float> smooth_facets_colors;
    // std::vector<float> smooth_vertices;
    // std::vector<float> smooth_vertices_colors;
    // std::vector<int>   smooth_facets;
    // std::vector<float> smooth_normals;
};

using namespace CGAL::Three;
class Scene_c3t3_color_item
    : public Scene_c3t3_item
{
    Q_OBJECT
public:
    typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

    Scene_c3t3_color_item(bool is_surface = false)
        : Scene_c3t3_item(is_surface)
    { }
    Scene_c3t3_color_item(const C3t3& c3t3, bool is_surface = false)
        : Scene_c3t3_item(c3t3, is_surface)
    { }
    virtual ~Scene_c3t3_color_item() { }

    const Color_infos& color_infos() const { return _color_infos; }
    Color_infos& color_infos() { return _color_infos; }

    bool supportsRenderingMode(RenderingMode m) const  override{
        return (m != PointsPlusNormals && m != Points && m != ShadedPoints);
    }

    void c3t3_changed();
    virtual bool do_take_edge(const T3::Edge& edge) const { return true; }
    void draw(CGAL::Three::Viewer_interface* viewer) const override;
    void drawEdges(CGAL::Three::Viewer_interface* viewer) const override;

    void invalidateOpenGLBuffers() override;

    void computeColorElements(bool flat) const;

    void setRenderingMode2(RenderingMode m) {
        Scene_c3t3_item::setRenderingMode(m);
    }

public Q_SLOTS:
    void computeElements() const override;

protected:

    Color_infos _color_infos; // indexed on triangulation's elements (vertices, edges, facets, cells)
    Draw_Color_infos _color_infos_draw; // indexed on vertices of triangulation's triangles to draw
};

QVector4D cgal_plane_to_vector4d(EPICK::Plane_3 plane) {
  return {
    static_cast<float>(-plane.a()),
    static_cast<float>(-plane.b()),
    static_cast<float>(-plane.c()),
    static_cast<float>(-plane.d()) };
}

void Scene_c3t3_color_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
    Scene_c3t3_color_item* ncthis = const_cast<Scene_c3t3_color_item*>(this);
    if (renderingMode() == Gouraud)
    {
        ncthis->setRenderingMode2(Flat);
        Scene_c3t3_item::drawEdges(viewer);
        ncthis->setRenderingMode2(Gouraud);
    }
    else if (renderingMode() == GouraudPlusEdges)
    {
        ncthis->setRenderingMode2(FlatPlusEdges);
        Scene_c3t3_item::drawEdges(viewer);
        ncthis->setRenderingMode2(GouraudPlusEdges);
    }
    else
    {
        Scene_c3t3_item::drawEdges(viewer);
    }
}

void Scene_c3t3_color_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
    Scene_c3t3_color_item* ncthis = const_cast<Scene_c3t3_color_item*>(this);
    bool flat_rendering_mode = (renderingMode() == Flat || renderingMode() == FlatPlusEdges);
    if (this->getBuffersFilled() && _color_infos_draw.flat_rendering_mode != flat_rendering_mode)
    {
        ncthis->_color_infos_draw.flat_rendering_mode = flat_rendering_mode;
        getTriangleContainer(T3_faces)->reset_vbos(COLORS);
        computeColorElements(flat_rendering_mode);
        this->initializeBuffers(viewer);
    }
    if (renderingMode() == Gouraud)
    {
        ncthis->setRenderingMode2(Flat);
        Scene_c3t3_item::draw(viewer);
        ncthis->setRenderingMode2(Gouraud);
    }
    else if (renderingMode() == GouraudPlusEdges)
    {
        ncthis->setRenderingMode2(FlatPlusEdges);
        Scene_c3t3_item::draw(viewer);
        ncthis->setRenderingMode2(GouraudPlusEdges);
    }
    else
    {
        Scene_c3t3_item::draw(viewer);
    }
}
void Scene_c3t3_color_item::invalidateOpenGLBuffers()
{
    Scene_c3t3_item::invalidateOpenGLBuffers();
}

void Scene_c3t3_color_item::computeColorElements(bool flat) const
{
    if (flat)
    {
        getTriangleContainer(T3_faces)->allocate(
            Triangle_container::FColors,
            const_cast<std::vector<float>&>(_color_infos_draw.flat_facets_colors).data(),
            static_cast<int>(_color_infos_draw.flat_facets_colors.size()*sizeof(float)));
    }
    else
    {
        getTriangleContainer(T3_faces)->allocate(
            Triangle_container::FColors,
            const_cast<std::vector<float>&>(_color_infos_draw.smooth_facets_colors).data(),
            static_cast<int>(_color_infos_draw.smooth_facets_colors.size()*sizeof(float)));
    }
}
void Scene_c3t3_color_item::computeElements() const
{
    Scene_c3t3_item::computeElements();
    Scene_c3t3_color_item* ncthis = const_cast<Scene_c3t3_color_item*>(this);
    ncthis->_color_infos_draw.flat_rendering_mode = (renderingMode() == Flat || renderingMode() == FlatPlusEdges);
    computeColorElements(ncthis->_color_infos_draw.flat_rendering_mode);
}

void Scene_c3t3_color_item::c3t3_changed()
{
    Scene_c3t3_item::c3t3_changed();

    // 1. build map
    std::map<C3t3::Vertex_handle, int> vertex_to_color_index;
    std::size_t vertices_index = 0;
    for (const C3t3::Vertex_handle& vertex : c3t3().triangulation().all_vertex_handles())
    {
        vertex_to_color_index[vertex] = vertices_index;
        ++vertices_index;
    }
    // 2. build colors and smooth infos
    _color_infos_draw.flat_facets_colors.clear();
    // _color_infos_draw.smooth_vertices_position.clear();
    _color_infos_draw.smooth_facets_colors.clear();
    std::size_t facets_index = -1;
    for (const C3t3::Facet& facet : c3t3().triangulation().all_facets())
    {
        ++facets_index;
        if (!this->do_take_facet(facet))
            continue;
        // flat color
        float color_r = _color_infos.facets_colors[3*facets_index];
        float color_g = _color_infos.facets_colors[3*facets_index+1];
        float color_b = _color_infos.facets_colors[3*facets_index+2];
        _color_infos_draw.flat_facets_colors.push_back(color_r);_color_infos_draw.flat_facets_colors.push_back(color_g);_color_infos_draw.flat_facets_colors.push_back(color_b);
        _color_infos_draw.flat_facets_colors.push_back(color_r);_color_infos_draw.flat_facets_colors.push_back(color_g);_color_infos_draw.flat_facets_colors.push_back(color_b);
        _color_infos_draw.flat_facets_colors.push_back(color_r);_color_infos_draw.flat_facets_colors.push_back(color_g);_color_infos_draw.flat_facets_colors.push_back(color_b);
        // smooth color
        std::array<int, 3> v_indices = {
            vertex_to_color_index[facet.first->vertex((facet.second+1)&3)],
            vertex_to_color_index[facet.first->vertex((facet.second+2)&3)],
            vertex_to_color_index[facet.first->vertex((facet.second+3)&3)]
        };
        if(this->is_facet_oriented(facet))
            std::swap(v_indices[0], v_indices[1]);
        for (int i = 0; i < 3; i++)
        {
            // _color_infos_draw.smooth_vertices_position.push_back();
            const int& v_index = v_indices[i];
            float color_r = _color_infos.vertices_colors[3*v_index];
            float color_g = _color_infos.vertices_colors[3*v_index+1];
            float color_b = _color_infos.vertices_colors[3*v_index+2];
            _color_infos_draw.smooth_facets_colors.push_back(color_r);_color_infos_draw.smooth_facets_colors.push_back(color_g);_color_infos_draw.smooth_facets_colors.push_back(color_b);
        }
    }
}

namespace Io_color_mesh
{

template <typename Tr>
struct file_color_triangulation_infos
{
    using Point_3             = typename Tr::Point;
    using Surface_patch_index = typename Tr::Cell::Surface_patch_index;
    using Subdomain_index     = typename Tr::Cell::Subdomain_index;
    using Color               = CGAL::IO::Color;
    using Edge                = std::array<int, 2>;
    using Facet               = std::array<int, 3>;
    using Cell                = std::array<int, 4>;

    std::vector<Point_3> points;
    std::vector<Edge> edges;
    std::vector<Facet> facets;
    std::vector<Cell> cells;

    std::vector<Surface_patch_index> facets_subdomains;
    std::vector<Subdomain_index>     cells_subdomains;

    std::vector<Color> vertices_colors;
    std::vector<Color> edges_colors;
    std::vector<Color> facets_colors;
    std::vector<Color> cells_colors;

    bool is_CGAL_mesh;
};

void order_facet_indices(std::array<int, 3>& facet)
{
    // find the circular permutation that puts the smallest index in the first place.
    // strategy : check all permutation
    int n0 = (std::min)({facet[0],facet[1], facet[2]});
    do
    {
        std::rotate(std::begin(facet), std::next(std::begin(facet)), std::end(facet));
    }
    while(facet[0] != n0);
}

void order_cell_indices(std::array<int, 4>& cell)
{
    // find the permutation that puts the smallest index in the first place, and the biggest index in the last place.
    // strategy : swap elements. The cell order is conserved if the number of swaps is even.
    //            the biggest element is swaped to the end and the smallest to the start.
    int min_id = std::numeric_limits<int>::max();
    int min_pos = 0;
    int max_id = std::numeric_limits<int>::min();
    int max_pos = 0;
    for (int v = 0; v < 4; v++)
    {
        const int& v_id = cell[v];
        if (v_id < min_id)
        {
            min_id = v_id;
            min_pos = v;
        }
        if (v_id < min_id)
        {
            min_id = v_id;
            min_pos = v;
        }
    }
    int swap_nb = 0;
    if (max_pos != 3)
    {
        std::swap(cell[max_pos], cell[3]);
        swap_nb++;
    }
    if (min_pos != 0)
    {
        std::swap(cell[min_pos], cell[0]);
        swap_nb++;
    }
    if (swap_nb == 1)
    {
        std::swap(cell[1], cell[2]);
    }
}

template <typename Tr>
bool get_triangulation_color_file_infos(std::istream& is,
                                        file_color_triangulation_infos<Tr>& file_infos)
{
    // Copied from CGAL::SMDS::build_triangulation_from_file
    using Point_3             = typename file_color_triangulation_infos<Tr>::Point_3;
    using Surface_patch_index = typename file_color_triangulation_infos<Tr>::Surface_patch_index;
    using Subdomain_index     = typename file_color_triangulation_infos<Tr>::Subdomain_index;
    using Color               = typename file_color_triangulation_infos<Tr>::Color;
    using Edge                = typename file_color_triangulation_infos<Tr>::Edge;
    using Facet               = typename file_color_triangulation_infos<Tr>::Facet;
    using Cell                = typename file_color_triangulation_infos<Tr>::Cell;

    if(!is)
        return false;

    int dim;
    std::string word;

    is >> word >> dim; // MeshVersionFormatted 1
    is >> word >> dim; // Dimension 3

    CGAL_assertion(dim == 3);

    file_infos.is_CGAL_mesh = false;

    std::string line;
    while(std::getline(is, line) && line != "End")
    {
        // remove trailing whitespace, in particular a possible '\r' from Windows
        // end-of-line encoding
        if(!line.empty() && std::isspace(line.back())) {
            line.pop_back();
        }
        if (line.size() > 0 && line.at(0) == '#' &&
            line.find("CGAL::Mesh_complex_3_in_triangulation_3") != std::string::npos)
        {
            file_infos.is_CGAL_mesh = true; // with CGAL meshes, domain 0 should be kept
            continue;
        }

        if(line == "Vertices")
        {
            int number_of_vertices, vertex_dimension_ignored;
            is >> number_of_vertices;
            file_infos.points.reserve(number_of_vertices);
            for(int i=0; i < number_of_vertices; ++i)
            {
                typename Tr::Geom_traits::FT x,y,z;
                if(!(is >> x >> y >> z >> vertex_dimension_ignored))
                    return false;
                file_infos.points.emplace_back(x,y,z);
            }
        }
        if(line == "Vertices_colors")
        {
            int number_of_vertices;
            is >> number_of_vertices;
            file_infos.vertices_colors.reserve(number_of_vertices);
            for(int i=0; i < number_of_vertices; ++i)
            {
                int r, g, b;
                if(!(is >> r >> g >> b))
                    return false;
                Color color(static_cast<unsigned char>(r), static_cast<unsigned char>(g), static_cast<unsigned char>(b));
                file_infos.vertices_colors.push_back(color);
            }
        }
        if(line == "Edges")
        {
            int number_of_edges, edge_line_index_ignored;
            is >> number_of_edges;
            file_infos.edges.reserve(number_of_edges);
            for(int i=0; i < number_of_edges; ++i)
            {
                int n[2];
                if(!(is >> n[0] >> n[1] >> edge_line_index_ignored))
                    return false;
                Edge edge;
                edge[0] = n[0] - 1;
                edge[1] = n[1] - 1;
                file_infos.edges.push_back(edge);
            }
        }
        if(line == "Edges_colors")
        {
            int number_of_edges;
            is >> number_of_edges;
            file_infos.edges_colors.reserve(number_of_edges);
            for(int i=0; i < number_of_edges; ++i)
            {
                int r, g, b;
                if(!(is >> r >> g >> b))
                    return false;
                Color color(static_cast<unsigned char>(r), static_cast<unsigned char>(g), static_cast<unsigned char>(b));
                file_infos.edges_colors.push_back(color);
            }
        }
        if(line == "Triangles")
        {
            bool has_negative_surface_patch_ids = false;
            Surface_patch_index max_surface_patch_id = 0;
            int number_of_facets;
            is >> number_of_facets;
            file_infos.facets.reserve(number_of_facets);
            for(int i=0; i<number_of_facets; ++i)
            {
                int n[3];
                Surface_patch_index surface_patch_id;
                if(!(is >> n[0] >> n[1] >> n[2] >> surface_patch_id))
                    return false;
                has_negative_surface_patch_ids |= (surface_patch_id < 0);
                max_surface_patch_id = (std::max)(max_surface_patch_id, surface_patch_id);
                Facet facet;
                facet[0] = n[0] - 1;
                facet[1] = n[1] - 1;
                facet[2] = n[2] - 1;
                // find the circular permutation that puts the smallest index in the first place.
                order_facet_indices(facet);
                file_infos.facets.push_back(facet);
                file_infos.facets_subdomains.push_back(surface_patch_id);
            }
            if(has_negative_surface_patch_ids)
            {
                for(Surface_patch_index& patch_id : file_infos.facets_subdomains) {
                    if(patch_id < 0)
                        patch_id = max_surface_patch_id - patch_id;
                }
            }
        }
        if(line == "Triangles_colors")
        {
            int number_of_facets;
            is >> number_of_facets;
            file_infos.facets_colors.reserve(number_of_facets);
            for(int i=0; i < number_of_facets; ++i)
            {
                int r, g, b;
                if(!(is >> r >> g >> b))
                    return false;
                Color color(static_cast<unsigned char>(r), static_cast<unsigned char>(g), static_cast<unsigned char>(b));
                file_infos.facets_colors.push_back(color);
            }
        }
        if(line == "Tetrahedra")
        {
            int number_of_cells;
            is >> number_of_cells;
            file_infos.cells.reserve(number_of_cells);
            file_infos.cells_subdomains.reserve(number_of_cells);
            for(int i=0; i < number_of_cells; ++i)
            {
                int n[4];
                Subdomain_index cell_domain;
                if(!(is >> n[0] >> n[1] >> n[2] >> n[3] >> cell_domain))
                    return false;
                Cell cell;
                cell[0] = n[0] - 1;
                cell[1] = n[1] - 1;
                cell[2] = n[2] - 1;
                cell[3] = n[3] - 1;
                file_infos.cells.push_back(cell);
                file_infos.cells_subdomains.push_back(cell_domain);
            }
        }
        if(line == "Tetrahedra_colors")
        {
            int number_of_cells;
            is >> number_of_cells;
            file_infos.cells_colors.reserve(number_of_cells);
            for(int i=0; i < number_of_cells; ++i)
            {
                int r, g, b;
                if(!(is >> r >> g >> b))
                    return false;
                Color color(static_cast<unsigned char>(r), static_cast<unsigned char>(g), static_cast<unsigned char>(b));
                file_infos.cells_colors.push_back(color);
            }
        }
    }
    return true;
}

template<class Tr>
bool build_triangulation_from_color_file(std::istream& is,
                                         Tr& tr,
                                         Color_infos& color_infos,
                                         const bool verbose,
                                         const bool replace_domain_0,
                                         const bool allow_non_manifold)
{
    // read file infos
    file_color_triangulation_infos<Tr> file_infos;
    if (!get_triangulation_color_file_infos(is, file_infos))
        return false;

    if (verbose)
    {
        std::cout << file_infos.points.size() << " points" << std::endl;
        std::cout << file_infos.edges.size()  << " edges" << std::endl;
        std::cout << file_infos.facets.size() << " facets" << std::endl;
        std::cout << file_infos.cells.size()  << " cells" << std::endl;
    }

    // check infos
    if (file_infos.cells.empty())
        return false;
    CGAL_assertion(file_infos.cells.size() == file_infos.cells_subdomains.size());
    CGAL_assertion(file_infos.facets.size() == file_infos.facets_subdomains.size());

    // remove subdomain 0 cells and their colors
    int cells_number_to_add = file_infos.cells.size() - file_infos.cells_colors.size();
    if (cells_number_to_add < 0)
    {
        file_infos.cells_colors.erase(file_infos.cells_colors.end() + cells_number_to_add, file_infos.cells_colors.end());
    }
    else
    {
        file_infos.cells_colors.resize(file_infos.cells_colors.size() + cells_number_to_add, CGAL::IO::Color(0,0,0));
    }
    std::vector<std::array<int, 4>>::iterator cell_it = file_infos.cells.begin();
    typename std::vector<typename Tr::Cell::Subdomain_index>::iterator subdomain_it = file_infos.cells_subdomains.begin();
    std::vector<CGAL::IO::Color>::iterator cell_color_it = file_infos.cells_colors.begin();
    for (; cell_it != file_infos.cells.end();)
    {
        if (*subdomain_it == 0)
        {
            cell_it = file_infos.cells.erase(cell_it);
            subdomain_it = file_infos.cells_subdomains.erase(subdomain_it);
            cell_color_it = file_infos.cells_colors.erase(cell_color_it);
        }
        else
        {
            ++cell_it;
            ++subdomain_it;
            ++cell_color_it;
        }
    }

    // build triangulation
    using Facet               = typename file_color_triangulation_infos<Tr>::Facet;
    using Surface_patch_index = typename file_color_triangulation_infos<Tr>::Surface_patch_index;
    using Vertex_handle       = typename Tr::Vertex_handle;
    boost::unordered_map<Facet, Surface_patch_index> border_facets;
    for (std::size_t f = 0, number_of_facets = file_infos.facets.size(); f < number_of_facets; f++)
    {
        border_facets.emplace(file_infos.facets[f], file_infos.facets_subdomains[f]);
    }
    std::vector<Vertex_handle> vertex_handle_vector; // contains vertex handles in the same order as points, except that the first one is the infinite vertex
    if (!CGAL::SMDS_3::build_triangulation_impl(tr, file_infos.points, file_infos.cells, file_infos.cells_subdomains, border_facets,
                                                vertex_handle_vector,
                                                verbose, replace_domain_0,
                                                allow_non_manifold))
        return false;

    // assign color
    using Color = typename file_color_triangulation_infos<Tr>::Color;
    //   assign vertices colors
    std::map<Vertex_handle, std::size_t> triangulation_vertex_to_id; // map triangulation vertex to the file's vertex index
    //   skip the first vertex handle : it is the infinite vertex and is not present in file
    for (std::size_t v = 1, number_of_vertices = vertex_handle_vector.size(); v < number_of_vertices; v++)
    {
        triangulation_vertex_to_id.emplace(vertex_handle_vector[v], v-1);
    }
    color_infos.vertices_colors.reserve(3 * vertex_handle_vector.size());
    typename std::map<Vertex_handle, std::size_t>::const_iterator find_vertex_iterator;
    for (const Vertex_handle& vh : tr.all_vertex_handles())
    {
        find_vertex_iterator = triangulation_vertex_to_id.find(vh);
        if (find_vertex_iterator == triangulation_vertex_to_id.end())
        {
            // color not found or infinite vertex
            color_infos.vertices_colors.push_back(0);
            color_infos.vertices_colors.push_back(0);
            color_infos.vertices_colors.push_back(0);
        }
        else
        {
            const Color& color = file_infos.vertices_colors[find_vertex_iterator->second];
            color_infos.vertices_colors.push_back(static_cast<float>(color.r()) / 255.0f);
            color_infos.vertices_colors.push_back(static_cast<float>(color.g()) / 255.0f);
            color_infos.vertices_colors.push_back(static_cast<float>(color.b()) / 255.0f);
        }
    }
    //   assign triangles colors
    using Triangulation_facet = typename Tr::Facet;
    using Facet = typename file_color_triangulation_infos<Tr>::Facet;
    std::map<Facet, std::size_t> facet_to_id; // map facet to the file's facet index
    for (std::size_t f = 0, number_of_facets = file_infos.facets.size(); f < number_of_facets; f++)
    {
        facet_to_id.emplace(file_infos.facets[f], f);
    }
    color_infos.facets_colors.reserve(3 * facet_to_id.size());
    typename std::map<Facet, std::size_t>::const_iterator find_facet_iterator;
    for (const Triangulation_facet tr_facet : tr.all_facets())
    {
        std::array<Vertex_handle, 3> tr_facet_vertices;
        tr_facet_vertices[0] = tr_facet.first->vertex((tr_facet.second+1)&3);
        tr_facet_vertices[1] = tr_facet.first->vertex((tr_facet.second+2)&3);
        tr_facet_vertices[2] = tr_facet.first->vertex((tr_facet.second+3)&3);
        Facet vertices_file_index;
        bool vertex_not_found = false;
        for (int v = 0; v < 3; v++)
        {
            find_vertex_iterator = triangulation_vertex_to_id.find(tr_facet_vertices[v]);
            if (find_vertex_iterator == triangulation_vertex_to_id.end())
            {
                // current triangle vertex is infinite or not in file
                vertex_not_found = true;
                break;
            }
            vertices_file_index[v] = find_vertex_iterator->second;
        }
        if (vertex_not_found)
        {
            color_infos.facets_colors.push_back(0);
            color_infos.facets_colors.push_back(0);
            color_infos.facets_colors.push_back(0);
            continue;
        }
        order_facet_indices(vertices_file_index);
        // find face index
        find_facet_iterator = facet_to_id.find(vertices_file_index);
        if (find_facet_iterator == facet_to_id.end())
        {
            // try find face index in other order
            std::swap(vertices_file_index[1], vertices_file_index[2]);
            find_facet_iterator = facet_to_id.find(vertices_file_index);
            if (find_facet_iterator == facet_to_id.end())
            {
                // face is not present in file despite all vertices being present, should probably never happen
                std::cerr << "warning : a face has been added while building triangulation" << std::endl;
                color_infos.facets_colors.push_back(0);
                color_infos.facets_colors.push_back(0);
                color_infos.facets_colors.push_back(0);
                continue;
            }
        }
        const Color& color = file_infos.facets_colors[find_facet_iterator->second];
        color_infos.facets_colors.push_back(static_cast<float>(color.r()) / 255.0f);
        color_infos.facets_colors.push_back(static_cast<float>(color.g()) / 255.0f);
        color_infos.facets_colors.push_back(static_cast<float>(color.b()) / 255.0f);
    }
    //   assign tetrahedra colors
    using Cell_handle = typename Tr::Cell_handle;
    using Cell = typename file_color_triangulation_infos<Tr>::Cell;
    std::map<Cell, std::size_t> cell_to_id; // map cell to the file's cell index. Cells indices are ordered.
    for (std::size_t c = 0, number_of_cells = file_infos.cells.size(); c < number_of_cells; c++)
    {
        Cell cell_copy = file_infos.cells[c];
        order_cell_indices(cell_copy);
        cell_to_id.emplace(cell_copy, c);
    }
    color_infos.cells_colors.reserve(3 * cell_to_id.size());
    typename std::map<Cell, std::size_t>::const_iterator find_cell_iterator;
    for (const Cell_handle ch : tr.all_cell_handles())
    {
        std::array<Vertex_handle, 4> tr_cell_vertices;
        tr_cell_vertices[0] = ch->vertex(0);
        tr_cell_vertices[1] = ch->vertex(1);
        tr_cell_vertices[2] = ch->vertex(2);
        tr_cell_vertices[3] = ch->vertex(3);
        Cell vertices_file_index;
        bool vertex_not_found = false;
        for (int v = 0; v < 4; v++)
        {
            find_vertex_iterator = triangulation_vertex_to_id.find(tr_cell_vertices[v]);
            if (find_vertex_iterator == triangulation_vertex_to_id.end())
            {
                // current cell vertex is infinite or not in file
                vertex_not_found = true;
                break;
            }
            vertices_file_index[v] = find_vertex_iterator->second;
        }
        if (vertex_not_found)
        {
            color_infos.cells_colors.push_back(0);
            color_infos.cells_colors.push_back(0);
            color_infos.cells_colors.push_back(0);
            continue;
        }
        order_cell_indices(vertices_file_index);
        // find cell index
        find_cell_iterator = cell_to_id.find(vertices_file_index);
        if (find_cell_iterator == cell_to_id.end())
        {
            // try find face index in other order
            std::swap(vertices_file_index[1], vertices_file_index[2]);
            find_cell_iterator = cell_to_id.find(vertices_file_index);
            if (find_cell_iterator == cell_to_id.end())
            {
                // cell is not present in file despite all vertices being present, should probably never happen
                std::cerr << "warning : a cell has been added while building triangulation" << std::endl;
                color_infos.cells_colors.push_back(0);
                color_infos.cells_colors.push_back(0);
                color_infos.cells_colors.push_back(0);
                continue;
            }
        }
        const Color& color = file_infos.cells_colors[find_cell_iterator->second];
        color_infos.cells_colors.push_back(static_cast<float>(color.r()) / 255.0f);
        color_infos.cells_colors.push_back(static_cast<float>(color.g()) / 255.0f);
        color_infos.cells_colors.push_back(static_cast<float>(color.b()) / 255.0f);
    }

    return true;
}

bool write_color_c3t3_to_color_file(std::ostream& os,
                                    const C3t3& c3t3,
                                    const Color_infos& color_infos)
{
    return false;
}

} // end namespace Io_color_mesh

class CGAL_Lab_c3t3_color_io_plugin :
                                      public QObject,
                                      public CGAL::Three::CGAL_Lab_io_plugin_interface
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

public:
    QString name() const override { return "C3t3_io_colors_plugin"; }
    QString nameFilters() const override { return "ascii (*.mesh)"; }
    QString saveNameFilters() const override{
        return "ascii (*.mesh)"; }
    QString loadNameFilters() const override {
        return "ascii (*.mesh)"; }
    bool canLoad(QFileInfo) const override;
    QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override;

    bool canSave(const CGAL::Three::Scene_item*) override;
    bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& ) override;
    bool isDefaultLoader(const Scene_item* item) const override{
        if(qobject_cast<const Scene_c3t3_color_item*>(item))
            return true;
        return false;
    }

private:
    void update_c3t3_and_colors(C3t3& c3t3, Color_infos& color_infos);
};


bool CGAL_Lab_c3t3_color_io_plugin::canLoad(QFileInfo fi) const {
    if(fi.suffix().toLower() != "mesh")
        return false;
    std::ifstream in(fi.filePath().toUtf8(),
                     std::ios_base::in|std::ios_base::binary);
    if(!in) {
        std::cerr << "Error! Cannot open file "
                  << (const char*)fi.filePath().toUtf8() << std::endl;
        return false;
    }
    // file header must contain comment "# Color"
    std::string line;
    std::istringstream iss;
    std::getline(in, line); // MeshVersionFormatted 1
    std::getline(in, line); // Dimension 3
    std::getline(in, line); // Vertices or #...
    iss.str(line);
    std::string keyword;
    while (iss >> keyword && keyword == "#")
    {
        if (iss >> keyword && keyword == "Color")
            return true;
        std::getline(in, line);
        iss.clear();
        iss.str(line);
    }
    std::cout << "end : " << keyword << std::endl;
    std::cerr << "File is not color mesh" << std::endl;
    return false;
}

QList<Scene_item*>
CGAL_Lab_c3t3_color_io_plugin::load(
    QFileInfo fileinfo, bool& ok, bool add_to_scene) {

    // Open file
    ok = true;
    std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::in);
    if(!in) {
        std::cerr << "Error! Cannot open file "
                  << (const char*)fileinfo.filePath().toUtf8() << std::endl;
        ok = false;
        return QList<Scene_item*>();
    }
    Scene_c3t3_color_item* item = new Scene_c3t3_color_item();

    if(fileinfo.size() == 0)
    {
        CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
        item->setName(fileinfo.completeBaseName());
        if(add_to_scene)
            CGAL::Three::Three::scene()->addItem(item);
        return QList<Scene_item*>()<< item;
    }
    if (fileinfo.suffix().toLower() == "mesh")
    {
        CGAL_assertion(!(!in));

        item->setName(fileinfo.completeBaseName());
        item->set_valid(false);

        if(Io_color_mesh::build_triangulation_from_color_file(in, item->c3t3().triangulation(), item->color_infos(),
                                                        /*verbose = */false, /*replace_subdomain_0 = */false, /*allow_non_manifold = */true))
        {
            update_c3t3_and_colors(item->c3t3(), item->color_infos());

            item->resetCutPlane();
            item->c3t3_changed();
            if(add_to_scene)
                CGAL::Three::Three::scene()->addItem(item);
            return QList<Scene_item*>()<<item;
        }
        else if(item->c3t3().triangulation().number_of_finite_cells() == 0)
        {
            QMessageBox::warning((QWidget*)NULL, tr("C3t3_io_color_plugin"),
                                 tr("No finite cell provided.\n"
                                    "Nothing to display."),
                                 QMessageBox::Ok);
        }
    }

    // if all loading failed...
    delete item;
    ok = false;
    return QList<Scene_item*>();
}

bool CGAL_Lab_c3t3_color_io_plugin::canSave(const CGAL::Three::Scene_item* item)
{
    // This plugin supports c3t3 items.
    return qobject_cast<const Scene_c3t3_color_item*>(item);
}

bool
    CGAL_Lab_c3t3_color_io_plugin::
    save(QFileInfo fileinfo, QList<Scene_item *> &items)
{
    Scene_item* item = items.front();
    const Scene_c3t3_color_item* c3t3_item = qobject_cast<const Scene_c3t3_color_item*>(item);
    if ( NULL == c3t3_item )
    {
        return false;
    }

    QString path = fileinfo.absoluteFilePath();

    if (fileinfo.suffix() == "mesh")
    {
        std::ofstream medit_file (qPrintable(path));
        Io_color_mesh::write_color_c3t3_to_color_file(medit_file, c3t3_item->c3t3(), c3t3_item->color_infos());
        items.pop_front();
        return true;
    }
    else
        return false;
}

void
    CGAL_Lab_c3t3_color_io_plugin::
    update_c3t3_and_colors(C3t3& c3t3, Color_infos& color_infos)
{
    using Cell_handle = C3t3::Triangulation::Cell_handle;
    // uses c3t3.add_to_complex(), wich does not change the underlying triangulation
    //      so the color_infos still matches with the triangulation
    c3t3.rescan_after_load_of_triangulation(); //fix counters for facets and cells
    for (Cell_handle cit : c3t3.triangulation().finite_cell_handles())
    {
        CGAL_assertion(cit->subdomain_index() >= 0);
        if (cit->subdomain_index() != C3t3::Triangulation::Cell::Subdomain_index())
            c3t3.add_to_complex(cit, cit->subdomain_index());

        for (int i = 0; i < 4; ++i)
        {
            if (cit->surface_patch_index(i) > 0)
                c3t3.add_to_complex(cit, i, cit->surface_patch_index(i));
        }
    }

    //if there is no facet in the complex, we add the border facets.
    if (c3t3.number_of_facets_in_complex() == 0)
    {
        for (C3t3::Facet fit : c3t3.triangulation().finite_facets())
        {
            Cell_handle c = fit.first;
            Cell_handle nc = c->neighbor(fit.second);

            // By definition, Subdomain_index() is supposed to be the id of the exterior
            if (c->subdomain_index() != C3t3::Triangulation::Cell::Subdomain_index() &&
                nc->subdomain_index() == C3t3::Triangulation::Cell::Subdomain_index())
            {
                // Color the border facet with the index of its cell
                c3t3.add_to_complex(c, fit.second, c->subdomain_index());
            }
        }
    }

    // make sure that the colors sizes is the same as triangulation sizes
    int vertices_number_to_add = c3t3.triangulation().number_of_vertices() - (color_infos.vertices_colors.size() / 3);
    if (vertices_number_to_add < 0)
    {
        std::cerr << "warning : vertices colors does not match : " << vertices_number_to_add << std::endl;
        color_infos.vertices_colors.erase(color_infos.vertices_colors.end() + 3*vertices_number_to_add, color_infos.vertices_colors.end());
    }
    else
    {
        color_infos.vertices_colors.resize(color_infos.vertices_colors.size() + 3*vertices_number_to_add, 0);
    }
    int edges_number_to_add = c3t3.triangulation().number_of_edges() - (color_infos.edges_colors.size() / 3);
    if (edges_number_to_add < 0)
    {
        std::cerr << "warning : edges colors does not match : " << edges_number_to_add << std::endl;
        color_infos.edges_colors.erase(color_infos.edges_colors.end() + edges_number_to_add, color_infos.edges_colors.end());
    }
    else
    {
        color_infos.edges_colors.resize(color_infos.edges_colors.size() + 3*edges_number_to_add, 0);
    }
    int facets_number_to_add = c3t3.triangulation().number_of_facets() - (color_infos.facets_colors.size() / 3);
    if (facets_number_to_add < 0)
    {
        std::cerr << "warning : facets colors does not match : " << facets_number_to_add << std::endl;
        color_infos.facets_colors.erase(color_infos.facets_colors.end() + facets_number_to_add, color_infos.facets_colors.end());
    }
    else
    {
        color_infos.facets_colors.resize(color_infos.facets_colors.size() + 3*facets_number_to_add, 0);
    }
    int cells_number_to_add = c3t3.triangulation().number_of_cells() - (color_infos.cells_colors.size() / 3);
    if (cells_number_to_add < 0)
    {
        std::cerr << "warning : cells colors does not match : " << cells_number_to_add << std::endl;
        color_infos.cells_colors.erase(color_infos.cells_colors.end() + cells_number_to_add, color_infos.cells_colors.end());
    }
    else
    {
        color_infos.cells_colors.resize(color_infos.cells_colors.size() + 3*cells_number_to_add, 0);
    }
}

#include <QtPlugin>
#include "C3t3_io_colors_plugin.moc"
