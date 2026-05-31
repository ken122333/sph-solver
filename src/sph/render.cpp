#include "render.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

GLuint renderer::compile_shader(GLenum type, const char* src)
{
    GLuint sh = glCreateShader(type);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);

    GLint ok = 0;
    glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);

    if (!ok)
    {
        GLint len = 0;
        glGetShaderiv(sh, GL_INFO_LOG_LENGTH, &len);

        std::vector<char> log((size_t)len + 1);
        glGetShaderInfoLog(sh, len, nullptr, log.data());

        std::fprintf(stderr, "Shader compile error:\n%s\n", log.data());

        glDeleteShader(sh);
        return 0;
    }

    return sh;
}

GLuint renderer::link_program(GLuint vs, GLuint fs)
{
    GLuint prog = glCreateProgram();

    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);

    GLint ok = 0;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);

    if (!ok)
    {
        GLint len = 0;
        glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &len);

        std::vector<char> log((size_t)len + 1);
        glGetProgramInfoLog(prog, len, nullptr, log.data());

        std::fprintf(stderr, "Program link error:\n%s\n", log.data());

        glDeleteProgram(prog);
        return 0;
    }

    return prog;
}

float renderer::to_ndc(float x, float a, float b)
{
    return ((x - a) / (b - a) * 2.0f - 1.0f) * 0.9;
}

std::vector<float> renderer::make_domain_outline_ndc(const domain& d)
{
    const float x0 = to_ndc(d.xmin, d.xmin, d.xmax);
    const float x1 = to_ndc(d.xmax, d.xmin, d.xmax);
    const float y0 = to_ndc(d.ymin, d.ymin, d.ymax);
    const float y1 = to_ndc(d.ymax, d.ymin, d.ymax);

    return {
        x0, y0,  x1, y0,
        x1, y0,  x1, y1,
        x1, y1,  x0, y1,
        x0, y1,  x0, y0
    };
}

bool renderer::init(const solver& solver, const domain& dom)
{
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const char* vs_points = R"GLSL(
        #version 330 core
        layout (location = 0) in vec2 aPos;
        layout (location = 1) in vec3 aColor;
        uniform float uPointSizePx;
        uniform vec2  uAspectScale;
        out vec3 vColor;

        void main() {
            vec2 p = aPos * uAspectScale;
            gl_Position = vec4(p, 0.0, 1.0);
            gl_PointSize = uPointSizePx;
            vColor = aColor;
        }
    )GLSL";

    const char* fs_points = R"GLSL(
        #version 330 core
        out vec4 FragColor;
        in vec3 vColor;

        void main() {
            vec2 centered = gl_PointCoord * 2.0 - 1.0;
            float radius2 = dot(centered, centered);

            if (radius2 > 1.0)
                discard;

            FragColor = vec4(vColor, 1.0);
        }
    )GLSL";

    const char* vs_lines = R"GLSL(
        #version 330 core
        layout (location = 0) in vec2 aPos;
        uniform vec2 uAspectScale;

        void main() {
            vec2 p = aPos * uAspectScale;
            gl_Position = vec4(p, 0.0, 1.0);
        }
    )GLSL";

    const char* fs_lines = R"GLSL(
        #version 330 core
        out vec4 FragColor;
        uniform vec3 uColor;

        void main() {
            FragColor = vec4(uColor, 1.0);
        }
    )GLSL";

    GLuint pvs = compile_shader(GL_VERTEX_SHADER, vs_points);
    GLuint pfs = compile_shader(GL_FRAGMENT_SHADER, fs_points);
    GLuint lvs = compile_shader(GL_VERTEX_SHADER, vs_lines);
    GLuint lfs = compile_shader(GL_FRAGMENT_SHADER, fs_lines);

    if (!pvs || !pfs || !lvs || !lfs)
        return false;

    prog_points = link_program(pvs, pfs);
    prog_lines = link_program(lvs, lfs);

    glDeleteShader(pvs);
    glDeleteShader(pfs);
    glDeleteShader(lvs);
    glDeleteShader(lfs);

    if (!prog_points || !prog_lines)
        return false;

    p_uPointSize = glGetUniformLocation(prog_points, "uPointSizePx");
    p_uAspect = glGetUniformLocation(prog_points, "uAspectScale");

    l_uAspect = glGetUniformLocation(prog_lines, "uAspectScale");
    l_uColor = glGetUniformLocation(prog_lines, "uColor");

    glGenVertexArrays(1, &vaoP);
    glGenBuffers(1, &vboP);

    glBindVertexArray(vaoP);
    glBindBuffer(GL_ARRAY_BUFFER, vboP);

    glBufferData(
        GL_ARRAY_BUFFER,
        (GLsizeiptr)(solver.particles.size() * 5 * sizeof(float)),
        nullptr,
        GL_DYNAMIC_DRAW
    );

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0,
        2,
        GL_FLOAT,
        GL_FALSE,
        5 * (GLsizei)sizeof(float),
        (void*)0
    );

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(
        1,
        3,
        GL_FLOAT,
        GL_FALSE,
        5 * (GLsizei)sizeof(float),
        (void*)(2 * sizeof(float))
    );

    glBindVertexArray(0);

    std::vector<float> outline = make_domain_outline_ndc(dom);

    glGenVertexArrays(1, &vaoL);
    glGenBuffers(1, &vboL);

    glBindVertexArray(vaoL);
    glBindBuffer(GL_ARRAY_BUFFER, vboL);

    glBufferData(
        GL_ARRAY_BUFFER,
        (GLsizeiptr)(outline.size() * sizeof(float)),
        outline.data(),
        GL_STATIC_DRAW
    );

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0,
        2,
        GL_FLOAT,
        GL_FALSE,
        2 * (GLsizei)sizeof(float),
        (void*)0
    );

    glBindVertexArray(0);

    upload_particles(solver, dom);

    return true;
}

void renderer::resize_particle_buffer(const solver& solver)
{
    glBindBuffer(GL_ARRAY_BUFFER, vboP);

    glBufferData(
        GL_ARRAY_BUFFER,
        (GLsizeiptr)(solver.particles.size() * 5 * sizeof(float)),
        nullptr,
        GL_DYNAMIC_DRAW
    );
}

void renderer::upload_particles(const solver& solver, const domain& dom)
{
    scratch.clear();
    scratch.reserve(solver.particles.size() * 5);

    for (const auto& p : solver.particles)
    {
        float speed = std::sqrt(p.v.x * p.v.x + p.v.y * p.v.y);
        float t = std::min(speed / 4.0f, 1.0f);

        scratch.push_back(to_ndc(p.x.x, dom.xmin, dom.xmax));
        scratch.push_back(to_ndc(p.x.y, dom.ymin, dom.ymax));
        scratch.push_back(0.1f + 0.9f * t);
        scratch.push_back(0.8f * (1.0f - std::abs(t - 0.5f) * 2.0f));
        scratch.push_back(1.0f - 0.8f * t);
    }

    glBindBuffer(GL_ARRAY_BUFFER, vboP);

    glBufferSubData(
        GL_ARRAY_BUFFER,
        0,
        (GLsizeiptr)(scratch.size() * sizeof(float)),
        scratch.data()
    );
}

void renderer::render(
    const solver& solver,
    int framebuffer_width,
    int framebuffer_height
)
{
    glViewport(0, 0, framebuffer_width, framebuffer_height);

    float aspect = framebuffer_height > 0
        ? (float)framebuffer_width / (float)framebuffer_height
        : 1.0f;

    float sx = 1.0f;
    float sy = 1.0f;

    if (aspect >= 1.0f)
        sx = 1.0f / aspect;
    else
        sy = aspect;

    glClearColor(0.05f, 0.07f, 0.10f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw particles
    glUseProgram(prog_points);
    glUniform1f(p_uPointSize, point_size_px);
    glUniform2f(p_uAspect, sx, sy);

    glBindVertexArray(vaoP);
    glDrawArrays(GL_POINTS, 0, (GLsizei)solver.particles.size());
    glBindVertexArray(0);

    // Draw domain outline
    glUseProgram(prog_lines);
    glUniform2f(l_uAspect, sx, sy);
    glUniform3f(l_uColor, 0.9f, 0.9f, 0.9f);

    glBindVertexArray(vaoL);
    glLineWidth(2.0f);
    glDrawArrays(GL_LINES, 0, 8);
    glBindVertexArray(0);
}

void renderer::cleanup()
{
    glDeleteBuffers(1, &vboP);
    glDeleteVertexArrays(1, &vaoP);

    glDeleteBuffers(1, &vboL);
    glDeleteVertexArrays(1, &vaoL);

    glDeleteProgram(prog_points);
    glDeleteProgram(prog_lines);

    vaoP = 0;
    vboP = 0;
    vaoL = 0;
    vboL = 0;

    prog_points = 0;
    prog_lines = 0;
}
