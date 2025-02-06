#pragma once

#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>

#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include <exception>
#include <stdexcept>

namespace polyhedralGravity {

    /**
     * Namespace containing various mesh reading functionalities
     */
    namespace MeshReader {

        /**
         * Returns a polyhedral source consisting of vertices and faces by reading mesh input files.
         * @param fileNames files specifying a polyhedron
         * @return polyhedral source consisting of vertices and faces
         */
        std::tuple<std::vector<Array3>, std::vector<IndexArray3>> getPolyhedralSource(const std::vector<std::string> &fileNames);

        /**
         * Reads elements from a .obj file (Wavefront OBJ file format)
         * This method only supports vertex (v) and faces (f) as input.
         * @param filename of the input source without suffix
         * @see https://de.wikipedia.org/wiki/Wavefront_OBJ
         */
        std::tuple<std::vector<Array3>, std::vector<IndexArray3>> readObj(const std::string &filename);

        /**
         * Reads elements from a file format supported by Tegen (.node/.face, .off, .ply, .stl, .mesh)
         * @param fileNames two files (.node/.face) or one file supported by Tetgen
         * @return Polyhedral Source consisting of vertices and faces
         */
        std::tuple<std::vector<Array3>, std::vector<IndexArray3>> readTetgenFormat(const std::vector<std::string> &fileNames);
    };

}
