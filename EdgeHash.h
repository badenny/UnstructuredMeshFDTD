//
// Created by Bud Denny on 6/22/22.
//

#ifndef GMSHTEST_EDGEHASH_H
#define GMSHTEST_EDGEHASH_H

#include <array>
#include <functional>

using size_t = std::size_t;
using SArray = std::array<size_t, 2>;

struct EdgeHash {
  std::size_t operator()(SArray const &s) const noexcept {
    const std::size_t h1 = std::hash<std::size_t>{}(s[0]);
    const std::size_t h2 = std::hash<std::size_t>{}(s[1]);
    return h1 ^ (h2 << 1);
  }
};

#endif //GMSHTEST_EDGEHASH_H
