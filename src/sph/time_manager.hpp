#pragma once
#include <algorithm>
#include <cmath>

struct time_manager
{
    double dt_fixed = 0.0005;
    int    max_substeps = 16;
    double max_frame_time = 0.1;
    bool   paused = false;
    double time_scale = 1.0;
    double accumulator = 0.0;

    void begin_frame(double frame_dt)
    {
        if (paused) return;
        frame_dt = std::min(frame_dt, max_frame_time);
        accumulator += frame_dt * time_scale;
    }

    bool step(double& out_dt)
    {
        if (paused) return false;
        if (accumulator < dt_fixed) return false;
        if (max_substeps-- <= 0) { accumulator = 0.0; return false; }
        out_dt = dt_fixed;
        accumulator -= dt_fixed;
        return true;
    }

    void reset_budget(int steps = 8) { max_substeps = steps; }
};