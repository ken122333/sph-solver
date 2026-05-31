#pragma once
#include <cstdint>
#include "math/vec2.hpp"


enum class particle_type : uint8_t
{
    fluid = 0,
    boundary = 1,
    ghost = 2
};

struct particle
{
    vec2 x;      // position
    vec2 v;      // velocity
    vec2 v_half; // velocity at the half time step
    vec2 a;      // acceleration                                                                                                                    

    float rho = 0.0f;   // density
    float p = 0.0f;   // pressure
    float m = 1.0f;   // mass

    particle_type type = particle_type::fluid;
};
