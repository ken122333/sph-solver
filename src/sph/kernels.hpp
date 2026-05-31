#pragma once
#include "math/vec2.hpp"
#include <cmath>

namespace kernels {
    inline float W_poly6(float r)
    {
        if (r >= h) return 0.0f;
        const float x = (h * h - r * r);
        const float c = 4.0f / (3.1415926f * std::pow(h, 8));
        return c * x * x * x;
    }
 
    inline vec2 grad_W_spiky(vec2 rij, float r)
    {
        if (r <= 1e-6f || r >= h) return { 0,0 };
        const float c = -30.0f / (3.1415926f * std::pow(h, 5));
        const float f = c * (h - r) * (h - r) / r;
        return rij * f;
    }
 


    inline float lap_W_visc(float r)
    {
        if (r >= h) return 0.0f;
        const float c = 40.0f / (3.1415926f * std::pow(h, 5));
        return c * (h - r);
    }

};
