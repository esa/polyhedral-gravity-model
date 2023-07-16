#pragma once

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <string>
#include "GravityModelData.h"

namespace polyhedralGravity {

    /* Forward declaration of Polyhedron */
    class Polyhedron;

    /**
     * Variant of possible polyhedron sources (composed of members, read from file), but not the polyhedron itself
     * @example Utilized in the Python interface which does not expose the Polyhedron class
     */
    using PolyhedralSource = std::variant<std::tuple<std::vector<Array3>, std::vector<IndexArray3>>, std::vector<std::string>>;

    /**
     * Variant of possible polyhedral sources (direct, composed of members, read from file), including the polyhedron
     * @example Utilized in the C++ implementation of the polyhedral-gravity model
     */
    using PolyhedronOrSource = std::variant<Polyhedron, std::tuple<std::vector<Array3>, std::vector<IndexArray3>>, std::vector<std::string>>;

    /**
     * Data structure containing the model data of one polyhedron. This includes nodes, edges (faces) and elements.
     * The index always starts with zero!
    */
    class Polyhedron {

        /**
         * A vector containing the vertices of the polyhedron.
         * Each node is an array of size three containing the xyz coordinates.
         */
        const std::vector<Array3> _vertices;

        /**
         * A vector containing the faces (triangles) of the polyhedron.
         * Each face is an array of size three containing the indices of the nodes forming the face.
         * Since every face consist of three nodes, every face consist of three segments. Each segment consists of
         * two nodes.
         * @example face consisting of {1, 2, 3} --> segments: {1, 2}, {2, 3}, {3, 1}
         */
        const std::vector<IndexArray3> _faces;


    public:

        /**
         * Generates an empty polyhedron.
         */
        Polyhedron()
                : _vertices{},
                  _faces{} {}

        /**
         * Generates a polyhedron from nodes and faces.
         * @param nodes - vector containing the nodes
         * @param faces - vector containing the triangle faces.
         *
         * @note ASSERTS PRE-CONDITION that the in the indexing in the faces vector starts with zero!
         * @throws runtime_error if no face contains the node zero indicating mathematical index
         */
        Polyhedron(std::vector<Array3> nodes, std::vector<IndexArray3> faces)
                : _vertices{std::move(nodes)},
                  _faces{std::move(faces)} {
            //Checks that the node with index zero is actually used
            if (_faces.end() == std::find_if(_faces.begin(), _faces.end(), [&](auto &face) {
                return face[0] == 0 || face[1] == 0 || face[2] == 0;
            })) {
                throw std::runtime_error("The node with index zero (0) was never used in any face! This is "
                                         "no valid polyhedron. Probable issue: Started numbering the vertices of "
                                         "the polyhedron at one (1).");
            }
        }

        /**
         * Creates a polyhedron from a tuple of nodes and faces.
         * @param data - tuple of nodes and faces
         *
         * @note ASSERTS PRE-CONDITION that the in the indexing in the faces vector starts with zero!
         * @throws runtime_error if no face contains the node zero indicating mathematical index
         */
        explicit Polyhedron(std::tuple<std::vector<Array3>, std::vector<IndexArray3>> data)
                : Polyhedron(std::get<std::vector<Array3>>(data), std::get<std::vector<IndexArray3>>(data)) {}


        /**
         * Default destructor
         */
        ~Polyhedron() = default;

        /**
         * Returns the vertices of this polyhedron
         * @return vector of cartesian coordinates
         */
        [[nodiscard]] const std::vector<Array3> &getVertices() const {
            return _vertices;
        }

        /**
         * Returns the vertex at a specific index
         * @param index - size_t
         * @return cartesian coordinates of the vertex at index
         */
        [[nodiscard]] const Array3 &getVertex(size_t index) const {
            return _vertices[index];
        }

        /**
         * The number of points (nodes) that make up the polyhedron.
         * @return a size_t
         */
        [[nodiscard]] size_t countVertices() const {
            return _vertices.size();
        }

        /**
         * Returns the number of faces (triangles) that make up the polyhedral.
         * @return a size_t
         */
        [[nodiscard]] size_t countFaces() const {
            return _faces.size();
        }

        /**
         * Returns the triangular faces of this polyhedron
         * @return vector of triangular faces, where each element size_t references a vertex in the vertices vector
         */
        [[nodiscard]] const std::vector<IndexArray3> &getFaces() const {
            return _faces;
        }

        /**
         * Returns the resolved face with its concrete cartesian coordinates at the given index.
         * @param index - size_t
         * @return triplet of vertices' cartesian coordinates
         */
        [[nodiscard]] Array3Triplet getFace(size_t index) const {
            return {_vertices[_faces[index][0]], _vertices[_faces[index][1]], _vertices[_faces[index][2]]};
        }

    };

}
