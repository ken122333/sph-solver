#pragma once
#include <vector>
#include "particle.hpp"
#include "spatial_hash_grid.h"

// SPH parameters
static constexpr float rho0 = 1000.0f;   // rest density
static constexpr float h = 0.03f;      // smoothing length (~2*dx)
static constexpr float c0 = 10.0f;      // speed of sound
static constexpr float gamma = 7.0f;
static constexpr float nu = 0.03f;       // viscosity
static const vec2 g = { 0.0f, -9.81f };  // gravity
static constexpr float k = rho0 * c0 * c0 / gamma;


struct domain
{
    float xmin = 0.0f, xmax = 1.0f;
    float ymin = 0.0f, ymax = 1.0f;
};

class solver
{
public:
    std::vector<particle> particles;
    std::vector<particle> ghost_particles;
    std::vector<particle> search_particles;

    void spawn_uniform_grid(const domain& dom, float dx, float margin, float rho0, float mass);
    void update(float dt);
    void solve_eos_pressure(particle& p);
    void solve_density(particle& p);
    void integrate(particle& p, float dt);
    void spawn_boundaries(const domain& dom, float dx, float mass);
    void spawn_block_top_center(const domain& dom, float dx, float block_width,float block_height,float top_margin, float rho0, float mass);
    void create_ghost_particle_boundary(const particle& p, const domain& dom, float h, std::vector<particle>& particles);
    particle create_ghost_particle(const particle& p, vec2 pos, vec2 vel);
    void solve_acceleration(particle& p);
    spatial_grid grid = spatial_grid(h);

};
