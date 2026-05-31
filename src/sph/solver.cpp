#include "solver.hpp"
#include "kernels.hpp"

#include <cmath>
#include <algorithm>

// Calculates pressure using the tait equation

void solver::solve_eos_pressure(particle& p)
{
    if (p.type != particle_type::fluid)
    {
        p.p = 0.0f;
        return;
    }

    float ratio = p.rho / rho0;
    p.p = k * (std::pow(ratio, gamma) - 1.0f);
    // Avoid negative pressure to prevent attractive pressure forces

    if (p.p < 0.0f)
        p.p = 0.0f;
}

void solver::spawn_boundaries(const domain& dom, float dx, float mass)
{
    int layers = 2;
    float spacing = dx/2;

    auto add_boundary = [&](float x, float y)
        {
            particle p;

            p.x = { x, y };
            p.v = { 0.0f, 0.0f };
            p.v_half = p.v;
            p.a = { 0.0f, 0.0f };

            p.rho = rho0;
            p.p = 0.0f;
            p.m = mass;

            p.type = particle_type::boundary;

            particles.push_back(p);
        };

    // Bottom + top
    for (float x = dom.xmin - layers * spacing;
        x <= dom.xmax + layers * spacing;
        x += spacing)
    {
        for (int l = 1; l <= layers; ++l)
        {
            add_boundary(x, dom.ymin - l * spacing);
            add_boundary(x, dom.ymax + l * spacing);
        }
    }

    // Left + right
    for (float y = dom.ymin - layers * spacing;
        y <= dom.ymax + layers * spacing;
        y += spacing)
    {
        for (int l = 1; l <= layers; ++l)
        {
            add_boundary(dom.xmin - l * spacing, y);
            add_boundary(dom.xmax + l * spacing, y);
        }
    }
}
// Computes acceleration on a particle
void solver::solve_acceleration(particle& pi) {

    // Gravity
    vec2 acc = g;

    // Get nearby grid cells around the particle
    std::vector<int> neighborHashes = grid.get_neighbor_hashes(pi.x);

    for (int hash : neighborHashes) {
        const std::vector<int>* neighbor_particles = grid.get_cell(hash);

        if (neighbor_particles) {
            for (int pj_index : *neighbor_particles) {
                
                const particle& pj = search_particles[pj_index];
                // Calculate distance between neighbor particle and target particle
                vec2 rij = pi.x - pj.x;
                float r2 = len2(rij);
                float r = std::sqrt(r2);
                if (r < 1e-6f)
                    continue;
                //Ignore particles further than the smoothing length
                if (r < h) {
                    vec2 gradW = kernels::grad_W_spiky(rij, r);
                    // Symmetric SPH pressure term
                    float pressure_term = (pi.p / (pi.rho * pi.rho) + pj.p / (pj.rho * pj.rho));
                    /*if (pj.type == particle_type::fluid)
                    {
                        pressure_term =
                            pi.p / (pi.rho * pi.rho) +
                            pj.p / (pj.rho * pj.rho);

                    }
                   
                    else if (pj.type == particle_type::boundary)
                    {
                        pressure_term =
                            pi.p / (pi.rho * pi.rho);

                    }*/
                    
                    
                    // Pressure contribution
                    acc = acc - (gradW * (pj.m * pressure_term));
                    //if (pj.type == particle_type::ghost)
                    //    continue;
                    float avg_mu = nu; // (pi.mu + pj.mu) / 2.0f for dynamic viscosity
                    float dot_product = dot(rij, gradW);
                    float eta2 = 0.01f * h * h;
                    float scalar_term = 2.0f * (pj.m / pj.rho) * (avg_mu) / (r*r + eta2);
                    vec2 vij = pj.v - pi.v;
                    acc = acc - vij * (scalar_term * dot_product);
                }
            }
        }
    }
    // Store final acceleration
    pi.a = acc;
    return;
};

// Computes density on a particle
void solver::solve_density(particle& pi) {
    pi.rho = 0.0f;

    // Find neighboring grid cells
    std::vector<int> neighborHashes = grid.get_neighbor_hashes(pi.x);

    for (int hash : neighborHashes) {
        const std::vector<int>* neighbor_particles = grid.get_cell(hash);

        if (neighbor_particles) {
            for (int pj_index : *neighbor_particles) {

                // Calculate distance between neighbor particle and target particle
                const particle& pj = search_particles[pj_index];
                vec2 rij = pi.x - pj.x;
                float r2 = len2(rij);
                float r = std::sqrt(r2);

                //Ignore particles further than the smoothing length
                if (r < h) {
                    pi.rho += pj.m * kernels::W_poly6(r);
                }
            }
        }
    }
    return;
}
bool near_left_wall(const particle& p, const domain& dom, float h)
{
    return (p.x.x - dom.xmin) < h;
}

bool near_right_wall(const particle& p, const domain& dom, float h)
{
    return (dom.xmax - p.x.x) < h;
}

bool near_bottom_wall(const particle& p, const domain& dom, float h)
{
    return (p.x.y - dom.ymin) < h;
}

bool near_top_wall(const particle& p, const domain& dom, float h)
{
    return (dom.ymax - p.x.y) < h;
}


particle solver::create_ghost_particle(const particle& p, vec2 pos, vec2 vel) {
    particle ghost = p;
    ghost.x = pos;
    ghost.v = vel;
    ghost.v_half = vel;
    ghost.a = vec2{ 0.0f, 0.0f };
    ghost.type = particle_type::ghost;
    return ghost;

}
void solver::create_ghost_particle_boundary(const particle& p, const domain& dom, float h, std::vector<particle>& particles)
{

    if (near_bottom_wall(p, dom, h))
        particles.push_back(create_ghost_particle(p, vec2(p.x.x, 2 * dom.ymin - p.x.y), vec2(p.v.x, -p.v.y)));
    else if (near_top_wall(p, dom, h))
        particles.push_back(create_ghost_particle(p, vec2(p.x.x, 2 * dom.ymax - p.x.y), vec2(p.v.x, -p.v.y)));
    if (near_left_wall(p, dom, h))
        particles.push_back(create_ghost_particle(p, vec2(2 * dom.xmin - p.x.x, p.x.y), vec2(-p.v.x, p.v.y)));
    else if (near_right_wall(p, dom, h))
        particles.push_back(create_ghost_particle(p, vec2(2 * dom.xmax - p.x.x, p.x.y), vec2(-p.v.x, p.v.y)));

   
}

void solver::integrate(particle& p, float dt) {

    p.v_half = p.v_half + p.a * dt;
    p.x = p.x + p.v_half * dt;
    p.v = p.v_half + p.a * (0.5f * dt);

};
void solver::update(float dt)
{
    domain dom;
    dom.xmin = 0.0f;
    dom.xmax = 1.0f;
    dom.ymin = 0.0f;
    dom.ymax = 1.0f;
    // Rebuild spatial hash grid
    grid.clear();
    ghost_particles.clear();
    search_particles.clear();




    for (size_t i = 0; i < particles.size(); ++i)
    {
        if (particles[i].type == particle_type::fluid)
            create_ghost_particle_boundary(particles[i], dom, h, ghost_particles);
    }
    search_particles.reserve(particles.size() + ghost_particles.size());
    for (const particle& p : particles)
        search_particles.push_back(p);

    for (const particle& g : ghost_particles)
        search_particles.push_back(g);

    for (int i = 0; i < search_particles.size(); ++i) {
        grid.insert(i, search_particles[i].x);
    }


    const size_t N = particles.size();
    if (N == 0) return;
    // Calculate density for all particles
    for (size_t i = 0; i < N; ++i)
    {
        if (particles[i].type == particle_type::fluid)
            solve_density(particles[i]);
        else
            particles[i].rho = rho0;
    }

    // Calculate pressure for all particles
    for (size_t i = 0; i < N; ++i)
    {
        solve_eos_pressure(particles[i]);
    }
    float rho_min = FLT_MAX;
    float rho_max = 0.0f;

    
    // Calculate acceleration for all particles
    for (size_t i = 0; i < N; ++i)
    {
        // Ignore boundary particles
        if (particles[i].type == particle_type::ghost) continue;
        solve_acceleration(particles[i]);
        rho_min = std::min(rho_min, particles[i].rho);
        rho_max = std::max(rho_max, particles[i].rho);
    }
   // printf("rho min %.1f rho max %.1f\n", rho_min, rho_max);



 
    // Integrate
    for (auto& p : particles)
    {
        // Ignore boundary particles
        if (p.type == particle_type::boundary) continue;
        solver::integrate(p, dt);
        

    }
}

void solver::spawn_block_top_center(
    const domain& dom,
    float dx,
    float block_width,
    float block_height,
    float top_margin,
    float rho0,
    float mass)
{
    particles.clear();

    if (dx <= 0.0f || block_width <= 0.0f || block_height <= 0.0f)
        return;

    // Center the block horizontally
    const float x0 = 0.40f * (dom.xmin + dom.xmax) - 0.5f * block_width;
    const float x1 = x0 + block_width;

    // Place the block near the top of the domain
    const float y1 = dom.ymax - top_margin;
    const float y0 = y1 - block_height;

    if (x1 <= x0 || y1 <= y0)
        return;

    const int nx = std::max<int>(1, (int)std::floor((x1 - x0) / dx) + 1);
    const int ny = std::max<int>(1, (int)std::floor((y1 - y0) / dx) + 1);

    particles.reserve((size_t)nx * (size_t)ny);

    for (int j = 0; j < ny; ++j)
    {
        const float y = y0 + j * dx;

        for (int i = 0; i < nx; ++i)
        {
            const float x = x0 + i * dx;

            particle p;

            p.x = vec2{ x, y };
            p.v = vec2{ 0.0f, 0.0f };
            p.v_half = p.v;
            p.a = vec2{ 0.0f, 0.0f };

            p.rho = rho0;
            p.p = 0.0f;
            p.m = mass;
            p.type = particle_type::fluid;

            particles.push_back(p);
        }
    }
};

void solver::spawn_uniform_grid(const domain& dom, float dx, float margin, float rho0, float mass)
{
    particles.clear();
    if (dx <= 0.0f) return;

    const float x0 = dom.xmin + margin;
    const float x1 = dom.xmax - margin;
    const float y0 = dom.ymin + margin;
    const float y1 = dom.ymax - margin;

    if (x1 <= x0 || y1 <= y0) return;

    const int nx = std::max<int>(1, (int)std::floor((x1 - x0) / dx) + 1);
    const int ny = std::max<int>(1, (int)std::floor((y1 - y0) / dx) + 1);

    particles.reserve((size_t)nx * (size_t)ny);

    for (int j = 0; j < ny; ++j)
    {
        const float y = y0 + j * dx;
        for (int i = 0; i < nx; ++i)
        {
            const float x = x0 + i * dx;

            particle p;
            p.x = vec2{ x, y };
            p.v = vec2{ 0.0f, 0.0f };
            p.v_half = p.v;
            p.a = vec2{ 0.0f, 0.0f };

            p.rho = rho0;
            p.p = 0.0f;
            p.m = mass;
            p.type = particle_type::fluid;

            particles.push_back(p);
        }
    }
}
