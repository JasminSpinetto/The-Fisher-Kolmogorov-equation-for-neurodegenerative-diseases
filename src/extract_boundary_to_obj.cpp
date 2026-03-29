#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <set>

using namespace dealii;

// --- NEW FUNCTION: Generates a curved C-shape using a Catmull-Rom spline ---
std::vector<Point<2>> generate_curved_hole_points(unsigned int points_per_segment = 15)
{
    // --- THE FIX: Adjusted control points for a wider, shorter, concave shape ---
    const std::vector<Point<2>> control_points = {
        // Point<2>(60, 78), Point<2>(75, 88), Point<2>(90, 78), // Top curve
        // Point<2>(98, 70), Point<2>(80, 68), Point<2>(62, 70)  // Bottom curve (more concave)
        Point<2>(42, 72), Point<2>(70, 82), Point<2>(98, 75), // Top curve
        Point<2>(95, 67), Point<2>(70, 72), Point<2>(50, 65)  // Bottom curve (more concave)
    };

    std::vector<Point<2>> spline_points;
    const unsigned int n = control_points.size();

    for (unsigned int i = 0; i < n; ++i) {
        const Point<2> &p0 = control_points[(i - 1 + n) % n];
        const Point<2> &p1 = control_points[i];
        const Point<2> &p2 = control_points[(i + 1) % n];
        const Point<2> &p3 = control_points[(i + 2) % n];

        for (unsigned int j = 0; j < points_per_segment; ++j) {
            float t = (float)j / points_per_segment;
            float t2 = t * t;
            float t3 = t2 * t;

            Point<2> point = 0.5f * ( (2.0f * p1) +
                                     (-p0 + p2) * t +
                                     (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                                     (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3 );
            spline_points.push_back(point);
        }
    }
    return spline_points;
}

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.msh> <output.obj>" << std::endl;
        return 1;
    }
    const std::string input_filename = argv[1];
    const std::string output_filename = argv[2];

    try {
        Triangulation<2> tria;
        GridIn<2> grid_in;
        grid_in.attach_triangulation(tria);
        std::ifstream input_file(input_filename);
        if (!input_file) return 1;
        grid_in.read_msh(input_file);

        // --- Use the proven logic to find and FILTER the boundary loops ---
        
        // 1. Build adjacency list of boundary edges
        std::map<unsigned int, std::vector<unsigned int>> adj;
        // ... (this part is the same as the last working version)
        std::map<std::pair<unsigned int, unsigned int>, unsigned int> edge_counts;
        for (const auto &cell : tria.active_cell_iterators()) {
            for (unsigned int v = 0; v < cell->n_vertices(); ++v) {
                unsigned int v0 = cell->vertex_index(v);
                unsigned int v1 = cell->vertex_index((v + 1) % cell->n_vertices());
                if (v0 > v1) std::swap(v0, v1);
                edge_counts[{v0, v1}]++;
            }
        }
        for (const auto &pair : edge_counts) {
            if (pair.second == 1) {
                adj[pair.first.first].push_back(pair.first.second);
                adj[pair.first.second].push_back(pair.first.first);
            }
        }
        
        // 2. Trace all distinct loops
        std::vector<std::vector<unsigned int>> all_loops;
        std::set<unsigned int> visited_nodes;
        for (auto const& [start_node, neighbors] : adj) {
            if (visited_nodes.find(start_node) == visited_nodes.end()) {
                std::vector<unsigned int> current_loop;
                unsigned int current_node = start_node;
                unsigned int prev_node = dealii::numbers::invalid_unsigned_int;
                while(visited_nodes.find(current_node) == visited_nodes.end()) {
                    current_loop.push_back(current_node);
                    visited_nodes.insert(current_node);
                    bool found_next = false;
                    for(const auto& next_candidate : adj[current_node]) {
                        if (next_candidate != prev_node) {
                            prev_node = current_node;
                            current_node = next_candidate;
                            found_next = true;
                            break;
                        }
                    }
                    if (!found_next) break;
                }
                if (current_loop.size() > 2)
                    all_loops.push_back(current_loop);
            }
        }
        std::cout << "Found " << all_loops.size() << " boundary loops." << std::endl;

        // 3. Find the largest loop and DISCARD the rest
        if (all_loops.empty()) throw std::runtime_error("No boundary loops found.");
        auto largest_loop_it = std::max_element(all_loops.begin(), all_loops.end(),
            [](const auto& a, const auto& b){ return a.size() < b.size(); });
        const std::vector<unsigned int>& outer_boundary_indices = *largest_loop_it;
        
        std::cout << "Identified outer boundary with " << outer_boundary_indices.size() << " vertices." << std::endl;
        if (all_loops.size() > 1)
            std::cout << "Discarded " << all_loops.size() - 1 << " smaller inner hole(s)." << std::endl;

        const std::vector<Point<2>>& all_vertices = tria.get_vertices();
        
        // 4. Generate the new curved hole
        std::vector<Point<2>> hole_vertices = generate_curved_hole_points();
        
        // 5. Write the final .obj file
        std::ofstream out(output_filename);
        out << "# Wavefront OBJ file generated by deal.II extractor\n";
        
        out << "# --- All Original Mesh Vertices ---\n";
        for (const auto& vertex : all_vertices) {
            out << "v 0 " << vertex[0] << " " << vertex[1] << "\n";
        }
        
        const unsigned int original_vertex_count = all_vertices.size();
        out << "\n# --- Vertices for the new Inner Curved Hole ---\n";
        for (const auto& vertex : hole_vertices) {
            out << "v 0 " << vertex[0] << " " << vertex[1] << "\n";
        }
        
        out << "\n# --- Line Segments for the Outer Boundary Loop ---\n";
        for (size_t i = 0; i < outer_boundary_indices.size(); ++i) {
            out << "l " << outer_boundary_indices[i] + 1 << " " 
                << outer_boundary_indices[(i + 1) % outer_boundary_indices.size()] + 1 << "\n";
        }
            
        out << "\n# --- Line Segments for the Inner Hole Loop ---\n";
        const unsigned int hole_vertex_count = hole_vertices.size();
        const unsigned int start_index_hole = original_vertex_count + 1;
        for (unsigned int i = 0; i < hole_vertex_count; ++i) {
            out << "l " << start_index_hole + i << " " << start_index_hole + ((i + 1) % hole_vertex_count) << "\n";
        }
        
        std::cout << "Successfully wrote .obj file to " << output_filename << std::endl;
    } catch (...) { return 1; }
    return 0;
}

