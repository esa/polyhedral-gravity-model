#include "TetgenAdapter.h"

namespace polyhedralGravity {
    TetgenAdapter::TetgenAdapter(const std::vector<std::string> &fileNames)
        : _tetgenio{}, _fileNames{fileNames}, _vertices{}, _faces{}
    {}

    PolyhedralSource TetgenAdapter::getPolyhedralSource() {
        //1. Step: Read in from files
        for (const auto &fileName: _fileNames) {
            size_t pos = fileName.find_last_of('.');
            std::string name = fileName.substr(0, pos);
            std::string suffix = fileName.substr(pos + 1);
            try {
                _suffixToOperation.at(suffix)(name);
            } catch (const std::out_of_range &e) {
                throw std::runtime_error("The TetGen Adapter does not support the requested file. " + fileName + " has an unknown suffix.");
            }
        }

        //2. Convert tetgenio to Polyhedron
        return std::make_tuple(_vertices, _faces);
    }

    void TetgenAdapter::readNode(const std::string &filename) {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}.node", filename);
        this->checkIntegrity(filename, 'v');
        try {
            _tetgenio.load_node(const_cast<char *>(filename.c_str()));
            this->addVertices();
        } catch (...) {
            throw std::runtime_error(DEFAULT_EXCEPTION_MSG);
        }
    }

    void TetgenAdapter::readFace(const std::string &filename) {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}.face", filename);
        this->checkIntegrity(filename, 'f');
        try {
            _tetgenio.load_face(const_cast<char *>(filename.c_str()));
            this->addFacesByTrifaces();
        } catch (...) {
            std::string message{DEFAULT_EXCEPTION_MSG};
            message.append(" A second possible issue could be a wrong file order, e.g. the .face file was read before "
                           "the .node file. In this case just reverse the parameters in the input file list.");
            throw std::runtime_error(message);
        }
    }

    void TetgenAdapter::readOff(const std::string &filename) {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}.off", filename);
        this->checkIntegrity(filename, 'a');
        try {
            _tetgenio.load_off(const_cast<char *>(filename.c_str()));
            this->addVerticesAndFacesByTriangulation();
        } catch (...) {
            throw std::runtime_error(DEFAULT_EXCEPTION_MSG);
        }
    }

    void TetgenAdapter::readPly(const std::string &filename) {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}.ply", filename);
        this->checkIntegrity(filename, 'a');
        try {
            _tetgenio.load_ply(const_cast<char *>(filename.c_str()));
            this->addVerticesAndFacesByTriangulation();
        } catch (...) {
            throw std::runtime_error(DEFAULT_EXCEPTION_MSG);
        }
    }

    void TetgenAdapter::readStl(const std::string &filename) {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}.stl", filename);
        this->checkIntegrity(filename, 'a');
        try {
            _tetgenio.load_stl(const_cast<char *>(filename.c_str()));
            this->addVerticesAndFacesByTriangulation();
        } catch (...) {
            throw std::runtime_error(DEFAULT_EXCEPTION_MSG);
        }
    }

    void TetgenAdapter::readMesh(const std::string &filename) {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}.mesh", filename);
        this->checkIntegrity(filename, 'a');
        try {
            _tetgenio.load_medit(const_cast<char *>(filename.c_str()), 0);
            this->addVerticesAndFacesByTriangulation();
        } catch (...) {
            throw std::runtime_error(DEFAULT_EXCEPTION_MSG);
        }
    }

    void TetgenAdapter::checkIntegrity(const std::string &filename, char what) const {
        if ((what == 'v' || what == 'a') && _tetgenio.numberofpoints != 0) {
            throw std::runtime_error(
                    "The Polyhedron already has well defined nodes! The information of " + filename + ".node is redundant!");
        } else if ((what == 'f' || what == 'a') && (_tetgenio.numberoftrifaces != 0 || _tetgenio.numberoffacets != 0)) {
            throw std::runtime_error(
                    "The Polyhedron already has well defined faces! The information of " + filename + ".node is redundant!");
        }
        POLYHEDRAL_GRAVITY_LOG_DEBUG("No duplicate information given. Integrity good!");
    }

    void TetgenAdapter::addVertices() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Converting tetgen's vertices to C++ Polyhedron");
        _vertices.clear();
        _vertices.reserve(_tetgenio.numberofpoints);
        for (size_t i = 0; i < _tetgenio.numberofpoints * 3; i += 3) {
            _vertices.push_back({_tetgenio.pointlist[i],
                                 _tetgenio.pointlist[i + 1],
                                 _tetgenio.pointlist[i + 2]});
        }
    }

    void TetgenAdapter::addFacesByTrifaces() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Converting tetgen's trifaces to C++ Polyhedron");
        _faces.clear();
        _faces.reserve(_tetgenio.numberoftrifaces);
        for (size_t i = 0; i < _tetgenio.numberoftrifaces * 3; i += 3) {
            _faces.push_back({static_cast<size_t>(_tetgenio.trifacelist[i]),
                              static_cast<size_t>(_tetgenio.trifacelist[i + 1]),
                              static_cast<size_t>(_tetgenio.trifacelist[i + 2])});
        }
    }

    void TetgenAdapter::addVerticesAndFacesByTriangulation() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Converting arbitrarily shaped facets of the tetgenio structure to triangles by using Tetgen's tetrahedralize method");
        tetgenbehavior tetgenbehavior;
        tetgenbehavior.zeroindex = 1;
        tetrahedralize(&tetgenbehavior, &_tetgenio, &_tetgenio);
        this->addVertices();
        this->addFacesByTrifaces();
    }

}// namespace polyhedralGravity
