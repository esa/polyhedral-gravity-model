#pragma once

#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>

#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/model/Definitions.h"
#include "polyhedralGravity/util/UtilityContainer.h"
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
         * @throws std::invalid_argument exception if the file type is unsupported by the implementation
         */
        PolyhedralSource getPolyhedralSource(const std::vector<std::string> &fileNames);

        /**
         * Reads elements from a .obj file (Wavefront OBJ file format)
         * This method only supports vertex (v) and faces (f) as input.
         *
         * This is also the file format of polyhedrons in some datasets, e.g.,in
         * https://pds.nasa.gov/ds-view/pds/viewDataset.jsp?dsid=EAR-A-5-DDR-RADARSHAPE-MODELS-V2.0
         * However, the suffix for these files is .tab
         * @param filename of the input source without suffix
         * @see Refer to https://de.wikipedia.org/wiki/Wavefront_OBJ for further help with the format
         */
        PolyhedralSource readObj(const std::string &filename);

        /**
         * Reads elements from a file format supported by Tegen (.node/.face, .off, .ply, .stl, .mesh)
         *
         * Delegates the processing to a TetgenAdapter.
         * @param fileNames two files (.node/.face) or one file supported by Tetgen
         * @return Polyhedral Source consisting of vertices and faces
         */
        PolyhedralSource readTetgenFormat(const std::vector<std::string> &fileNames);
    };

}
