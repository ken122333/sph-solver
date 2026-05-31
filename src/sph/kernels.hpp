#pragma once
#include "math/vec2.hpp"
#include <cmath>

namespace kernels {
    inline constexpr float pi = 3.1415926f;
    inline constexpr float h2 = h * h;
    inline constexpr float h5 = h2 * h2 * h;
    inline constexpr float h8 = h2 * h2 * h2 * h2;
    inline constexpr float poly6_c = 4.0f / (pi * h8);
    inline constexpr float spiky_grad_c = -30.0f / (pi * h5);
    inline constexpr float visc_lap_c = 40.0f / (pi * h5);

    inline float W_poly6(float r)
    {
        if (r >= h) return 0.0f;
        const float x = h2 - r * r;
        return poly6_c * x * x * x;
    }

    inline float W_poly6_r2(float r2)
    {
        if (r2 >= h2) return 0.0f;
        const float x = h2 - r2;
        return poly6_c * x * x * x;
    }
 
    inline vec2 grad_W_spiky(vec2 rij, float r)
    {
        if (r <= 1e-6f || r >= h) return { 0,0 };
        const float f = spiky_grad_c * (h - r) * (h - r) / r;
        return rij * f;
    }
 


    inline float lap_W_visc(float r)
    {
        if (r >= h) return 0.0f;
        return visc_lap_c * (h - r);
    }

};
