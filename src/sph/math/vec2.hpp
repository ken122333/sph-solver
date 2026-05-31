#pragma once

struct vec2
{
    float x = 0.0f;
    float y = 0.0f;

    vec2() = default;
    vec2(float x_, float y_) : x(x_), y(y_) {}
};

inline vec2 operator+(vec2 a, vec2 b) { return { a.x + b.x, a.y + b.y }; }
inline vec2 operator-(vec2 a, vec2 b) { return { a.x - b.x, a.y - b.y }; }
inline vec2 operator*(vec2 a, float s) { return { a.x * s, a.y * s }; }
inline vec2 operator/(vec2 a, float s) { return { a.x / s, a.y / s }; }
inline vec2 operator-(vec2 a) { return {-a.x, -a.y}; }
inline float dot(vec2 a, vec2 b) { return a.x * b.x + a.y * b.y; }
inline float len2(vec2 a) { return dot(a, a); }
