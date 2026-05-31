#pragma once

#include <vector>
#include <glad/gl.h>

#include "sph/solver.hpp"
#include "sph/particle.hpp"

class renderer {
public:
    bool init(const solver& solver, const domain& dom);
    void upload_particles(const solver& solver, const domain& dom);
    void render(const solver& solver, int framebuffer_width, int framebuffer_height);
    void resize_particle_buffer(const solver& solver);
    void cleanup();

private:
    GLuint vaoP = 0;
    GLuint vboP = 0;

    GLuint vaoL = 0;
    GLuint vboL = 0;

    GLuint prog_points = 0;
    GLuint prog_lines = 0;

    GLint p_uPointSize = -1;
    GLint p_uAspect = -1;
    GLint p_uColor = -1;

    GLint l_uAspect = -1;
    GLint l_uColor = -1;

    float point_size_px = 12.0f;

    std::vector<float> scratch;

private:
    GLuint compile_shader(GLenum type, const char* src);
    GLuint link_program(GLuint vs, GLuint fs);

    float to_ndc(float x, float a, float b);
    std::vector<float> make_domain_outline_ndc(const domain& d);
};
