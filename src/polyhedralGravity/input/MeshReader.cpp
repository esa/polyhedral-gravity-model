#include "MeshReader.h"

#include "TetgenAdapter.h"
#include "polyhedralGravity/util/UtilityString.h"

namespace polyhedralGravity {

    namespace {
        std::string malformedLineMessage(const std::string &filename, size_t lineNumber, const std::string &line,
                                         const std::string &reason) {
            std::stringstream sstream{};
            sstream << "Malformed line " << lineNumber << " in file '" << filename << "': \"" << line << "\". " << reason;
            return sstream.str();
        }

        /**
         * Parses a single vertex index of a face definition. Accepts the plain index as well as the
         * v/vt/vn notation of the OBJ format, of which only the vertex part is relevant here.
         */
        size_t parseObjVertexIndex(const std::string &token, const std::string &filename, size_t lineNumber,
                                   const std::string &line) {
            const std::string vertexPart = token.substr(0, token.find('/'));
            size_t charactersParsed = 0;
            long long index = 0;
            try {
                index = std::stoll(vertexPart, &charactersParsed);
            } catch (const std::logic_error &) {
                throw std::runtime_error(malformedLineMessage(filename, lineNumber, line,
                                                              "The vertex index '" + token + "' is not an integer."));
            }
            if (charactersParsed != vertexPart.size()) {
                throw std::runtime_error(malformedLineMessage(filename, lineNumber, line,
                                                              "The vertex index '" + token + "' is not an integer."));
            }
            // Relative (negative) indices of the OBJ format are unsupported and must not wrap around when cast
            if (index < 0) {
                throw std::runtime_error(malformedLineMessage(filename, lineNumber, line,
                                                              "The vertex index '" + token + "' is negative. Relative vertex "
                                                              "indices are not supported."));
            }
            return static_cast<size_t>(index);
        }
    }

    PolyhedralSource MeshReader::getPolyhedralSource(const std::vector<std::string> &fileNames) {
        // Input Sanity Check if the files exists
        for (const auto &fileName: fileNames) {
            if (!std::filesystem::exists(fileName)) {
                throw std::runtime_error("File '" + fileName + "' does not exist.");
            }
        }
        switch (fileNames.size()) {
            case 0:
                throw std::runtime_error("No mesh file given");
            case 1:
                if (util::ends_with(fileNames[0], ".obj", ".tab")) {
                    return readObj(fileNames[0]);
                } else {
                    // The TetGen Adapter complains if the suffix is unknown
                    // No need to check for suffices here, as this happens later anyway
                    return readTetgenFormat(fileNames);
                }
            case 2:
                // The TetGen Adapter complains if the suffix is unknown
                // No need to check for suffices here, as this happens later anyway
                return readTetgenFormat(fileNames);
            default:
                throw std::invalid_argument("More than two mesh files given. There is no known mesh-format consisting of three files. "
                                            "The polyhedron will be over-specified!");
        }
    }

    PolyhedralSource MeshReader::readTetgenFormat(const std::vector<std::string> &fileNames) {
        return TetgenAdapter{fileNames}.getPolyhedralSource();
    }

    PolyhedralSource MeshReader::readObj(const std::string &filename) {
        std::vector<Array3> vertices{};
        std::vector<IndexArray3> faces{};
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the file {}", filename);
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file " + filename + " for reading.");
        }
        std::string line;
        size_t lineNumber = 0;
        while (std::getline(file, line)) {
            ++lineNumber;
            std::stringstream ss(line);
            std::string keyword;
            // Skips empty lines, whitespace-only lines and comments
            if (!(ss >> keyword) || keyword[0] == '#') {
                continue;
            }

            if (keyword == "v") {
                Array3 vertex{};
                if (!(ss >> vertex[0] >> vertex[1] >> vertex[2])) {
                    throw std::runtime_error(malformedLineMessage(filename, lineNumber, line,
                                                                  "A vertex (v) requires three numeric coordinates."));
                }
                vertices.push_back(vertex);
            } else if (keyword == "f") {
                IndexArray3 face{};
                std::string token;
                for (size_t i = 0; i < 3; ++i) {
                    if (!(ss >> token)) {
                        throw std::runtime_error(malformedLineMessage(filename, lineNumber, line,
                                                                      "A face (f) requires three vertex indices since only triangular meshes are supported."));
                    }
                    face[i] = parseObjVertexIndex(token, filename, lineNumber, line);
                }
                // A fourth index would silently be dropped and yield a different mesh than the file describes
                if (ss >> token) {
                    throw std::runtime_error(malformedLineMessage(filename, lineNumber, line,
                                                                  "A face (f) references more than three vertices, but only triangular meshes are supported."));
                }
                faces.push_back(face);
            }
            // Every other keyword (vn, vt, vp, g, o, s, usemtl, ...) carries no information for the gravity model
        }
        file.close();
        return {vertices, faces};
    }
}