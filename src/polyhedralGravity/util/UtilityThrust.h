#pragma once

#include <utility>
#include "thrust/iterator/zip_iterator.h"

namespace polyhedralGravity::util {

    /**
     * Returns a zip iterator for a pack of Iterators
     * @tparam Iterator - should be a iterator
     * @param args - can be any number of iterators
     * @return a thrust::zip_iterator
     */
    template<typename ...Iterator>
    auto zip(const Iterator&... args) {
        return thrust::make_zip_iterator(thrust::make_tuple(args...));
    }

    /**
     * Returns a zip iterator for a pack of Containers
     * @tparam Container - should be a container on which one can call std::begin()
     * @param args - can be any number of containers
     * @return a thrust::zip_iterator for the beginning of all containers
     */
    template<typename ...Container>
    auto zipBegin(const Container&... args) {
        return zip(std::begin(args)...);
    }

    /**
     * Returns a zip iterator for a pack of Containers
     * @tparam Container - should be a container on which one can call std::end()
     * @param args - can be any number of containers
     * @return a thrust::zip_iterator for the end of all containers
     */
    template<typename ...Container>
    auto zipEnd(const Container&... args) {
        return zip(std::end(args)...);
    }

    /**
     * Returns a zip iterator pair for a pack of Containers
     * @tparam Container - should be a container on which one can call std::begin() and std::end()
     * @param args - can be any number of containers
     * @return a pair of thrust::zip_iterator for the beginning and end of all containers
     */
    template<typename ...Container>
    auto zipPair(const Container&... args) {
        auto begin = zipBegin(args...);
        auto end = zipEnd(args...);
        return std::make_pair(begin, end);
    }

}