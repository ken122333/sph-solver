#include "spatial_hash_grid.h"
#include <cmath>

int spatial_grid::compute_hash(vec2 pos) const {
    // Convert particle coordinate to integer coordinate
    int ix = static_cast<int>(std::floor(pos.x / h));
    int iy = static_cast<int>(std::floor(pos.y / h));

    // Hashing with large primes
    return (ix * 83492791) ^ (iy * 2654435761);
}

// Adds particle index to grid cell containing its position
void spatial_grid::insert(int particle_index, vec2 pos) {
    data[compute_hash(pos)].push_back(particle_index);
}

// Removes all particles from the grid. Called before rebuilding the grid after every time step
void spatial_grid::clear() {
    data.clear();
}

void spatial_grid::reserve(size_t cell_count) {
    data.reserve(cell_count);
}

// Returns the list of particle indices stored in a grid cell
const std::vector<int>* spatial_grid::get_cell(int hash) const {
    auto it = data.find(hash);
    if (it != data.end()) {
        return &(it->second);
    }
    return nullptr;
}

// Returns hashes for the current cell and all surrounding cells
std::array<int, 9> spatial_grid::get_neighbor_hashes(vec2 pos) const {
    std::array<int, 9> hashes{};
    int ix = static_cast<int>(std::floor(pos.x / h));
    int iy = static_cast<int>(std::floor(pos.y / h));

    int index = 0;

    for (int x = ix - 1; x <= ix + 1; ++x) {
        for (int y = iy - 1; y <= iy + 1; ++y) {
            hashes[index++] = (x * 83492791) ^ (y * 2654435761);
        }
    }

    return hashes;
}
