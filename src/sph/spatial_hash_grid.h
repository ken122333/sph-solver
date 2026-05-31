#pragma once
#include <vector>
#include <unordered_map>
#include "math/vec2.hpp"

class spatial_grid {
public:
    spatial_grid(float cell_size) : h(cell_size) {}
    void insert(int particleIndex, vec2 pos);
    void clear();
    std::vector<int> get_neighbor_hashes(vec2 pos);

    const std::vector<int>* get_cell(int hash) const;

    int compute_hash(vec2 pos) const;

private:
    float h;
    std::unordered_map<int, std::vector<int>> data;
};