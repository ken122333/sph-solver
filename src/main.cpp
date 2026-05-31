#include <cstdio>
#include <algorithm>
#include <cmath>
#include <string>

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include "sph/solver.hpp"
#include "sph/particle.hpp"
#include "sph/time_manager.hpp"

#include "sph/render.hpp"

static void glfw_error_callback(int error, const char* description)
{
    std::fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

static bool key_pressed(GLFWwindow* w, int key)
{
    static bool prev[512] = {};

    bool down = glfwGetKey(w, key) == GLFW_PRESS;
    bool pressed = down && !prev[key];

    prev[key] = down;

    return pressed;
}

static bool cursor_to_domain(GLFWwindow* window, const domain& dom, vec2& out_pos)
{
    double cursor_x = 0.0;
    double cursor_y = 0.0;
    glfwGetCursorPos(window, &cursor_x, &cursor_y);

    int window_width = 0;
    int window_height = 0;
    glfwGetWindowSize(window, &window_width, &window_height);

    if (window_width <= 0 || window_height <= 0)
        return false;

    float ndc_x = (float)(cursor_x / (double)window_width) * 2.0f - 1.0f;
    float ndc_y = 1.0f - (float)(cursor_y / (double)window_height) * 2.0f;

    int framebuffer_width = 0;
    int framebuffer_height = 0;
    glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);

    float aspect = framebuffer_height > 0
        ? (float)framebuffer_width / (float)framebuffer_height
        : 1.0f;

    float sx = 1.0f;
    float sy = 1.0f;

    if (aspect >= 1.0f)
        sx = 1.0f / aspect;
    else
        sy = aspect;

    float sim_x = (ndc_x / sx) / 0.9f;
    float sim_y = (ndc_y / sy) / 0.9f;

    out_pos.x = dom.xmin + ((sim_x + 1.0f) * 0.5f) * (dom.xmax - dom.xmin);
    out_pos.y = dom.ymin + ((sim_y + 1.0f) * 0.5f) * (dom.ymax - dom.ymin);

    return true;
}

static void apply_mouse_interaction(GLFWwindow* window, solver& solver, const domain& dom, float dt)
{
    bool repel = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    bool attract = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

    if (!repel && !attract)
        return;

    vec2 cursor;
    if (!cursor_to_domain(window, dom, cursor))
        return;

    constexpr float radius = 0.3f;
    constexpr float radius2 = radius * radius;
    constexpr float repel_strength = 100.0f;
    constexpr float attract_strength = 100.0f;
    constexpr float grab_damping = 5.0f;

    for (auto& p : solver.particles)
    {
        if (p.type != particle_type::fluid)
            continue;

        vec2 to_particle = p.x - cursor;
        float d2 = len2(to_particle);

        if (d2 >= radius2)
            continue;

        float distance = std::sqrt(std::max(d2, 1.0e-8f));
        float falloff = 1.0f - distance / radius;
        falloff *= falloff;

        vec2 direction = to_particle / distance;
        vec2 delta_v = { 0.0f, 0.0f };

        if (repel)
            delta_v = delta_v + direction * (repel_strength * falloff * dt);

        if (attract)
        {
            vec2 to_cursor = -direction;
            delta_v = delta_v + to_cursor * (attract_strength * falloff * dt);
            delta_v = delta_v - p.v * (grab_damping * falloff * dt);
        }

        p.v = p.v + delta_v;
        p.v_half = p.v_half + delta_v;
    }
}
int particle_count = 0;
int main()
{
    glfwSetErrorCallback(glfw_error_callback);

    if (!glfwInit())
        return 1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window =
        glfwCreateWindow(1280, 800, "SPH Solver", nullptr, nullptr);

    if (!window)
    {
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (gladLoadGL(glfwGetProcAddress) == 0)
    {
        std::fprintf(stderr, "Failed to initialize GLAD\n");

        glfwDestroyWindow(window);
        glfwTerminate();

        return 1;
    }

    solver solver;

    domain dom;
    dom.xmin = 0.0f;
    dom.xmax = 1.0f;
    dom.ymin = 0.0f;
    dom.ymax = 1.0f;

    const float dx = h * 0.5;//0.025f;
    const float rho0 = 1000.0f;
    const float mass = rho0 * dx * dx;
    const float margin = 1.2f * dx;

    auto respawn = [&]()
        {
            
            solver.spawn_block_top_center(dom,
                dx,   // particle spacing
                0.5f,    // block width
                0.5f,    // block height
                0.05f,   // distance from top boundary
                1000.0f, // rest density
                mass     // particle mass
            );
           // solver.spawn_uniform_grid(dom, dx, margin, rho0, mass);
           // solver.spawn_boundaries(dom, dx,mass);

            std::printf("Spawned %zu particles\n", solver.particles.size());
        };

    respawn();

    renderer renderer;

    if (!renderer.init(solver, dom))
    {
        std::fprintf(stderr, "Failed to initialize renderer\n");

        glfwDestroyWindow(window);
        glfwTerminate();

        return 1;
    }

    time_manager stepper;
    stepper.dt_fixed = 0.0005f;
    double last_time = glfwGetTime();
    double fps_time = glfwGetTime();
    double print_timer = 0.0;

    int frame_count = 0;

    while (!glfwWindowShouldClose(window))
    {
        double current_time = glfwGetTime();

        frame_count++;

        if (current_time - fps_time >= 1.0)
        {
            double fps = frame_count / (current_time - fps_time);

            std::string title =
                "SPH Simulation - FPS: " + std::to_string((int)fps) +" - Particle Count: " +std::to_string((int)solver.particles.size());

            glfwSetWindowTitle(window, title.c_str());

            frame_count = 0;
            fps_time = current_time;
        }

        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GLFW_TRUE);

        if (key_pressed(window, GLFW_KEY_P))
            stepper.paused = !stepper.paused;

        if (key_pressed(window, GLFW_KEY_R))
        {
            respawn();
            renderer.resize_particle_buffer(solver);
        }

        double now = glfwGetTime();
        double frame_dt = now - last_time;
        last_time = now;

        stepper.reset_budget(16);
        stepper.begin_frame(frame_dt);

        double dt;

        while (stepper.step(dt))
        {
            apply_mouse_interaction(window, solver, dom, (float)dt);
            solver.update((float)dt);
        }

        renderer.upload_particles(solver, dom);

        print_timer += frame_dt;

        if (print_timer >= 1.0 && !solver.particles.empty())
        {
            print_timer = 0.0;

            const auto& p0 = solver.particles[0];
            
            std::printf(
                "p0: x=(%.3f, %.3f) v=(%.3f, %.3f) rho=%.1f p=%.1f\n",
                p0.x.x,
                p0.x.y,
                p0.v.x,
                p0.v.y,
                p0.rho,
                p0.p
            );
        }

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);

        renderer.render(solver, w, h);

        glfwSwapBuffers(window);
    }

    renderer.cleanup();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
