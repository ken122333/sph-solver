#pragma once
#include <array>
#include <vector>
#include <unordered_map>
#include "math/vec2.hpp"

class spatial_grid {
public:
    spatial_grid(float cell_size) : h(cell_size) {}
    void insert(int particleIndex, vec2 pos);
    void clear();
    void reserve(size_t cell_count);
    std::array<int, 9> get_neighbor_hashes(vec2 pos) const;

    const std::vector<int>* get_cell(int hash) const;

    int compute_hash(vec2 pos) const;

private:
    float h;
    std::unordered_map<int, std::vector<int>> data;
};
