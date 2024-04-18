#include "rayt.h"

bool hit_sphere(const vec3& center, float radius, const rayt::Ray& r) {
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = 2.0f * dot(r.direction(), oc);
    float c = dot(oc, oc) - pow2(radius);
    float D = b * b - 4 * a * c;
    return (D > 0);
}

vec3 color(const rayt::Ray& r) {
    if (hit_sphere(vec3(0, 0, -1), 0.5f, r)) {
        return vec3(0.0f, 1.0f, 0.0f);
    }

    vec3 d = normalize(r.direction());
    float t = 0.5f * (r.direction().getY() + 1.0f);
    return lerp(t, vec3(0.5f, 0.7f, 1.0f), vec3(1));
}

int main()
{
    int nx = 200;
    int ny = 100;
    std::unique_ptr<rayt::Image> image(std::make_unique<rayt::Image>(nx, ny));

    vec3 x(4.0f, 0.0f, 0.0f);
    vec3 y(0.0f, 2.0f, 0.0f);
    vec3 z(-2.0f, -1.0f, -1.0f);
    std::unique_ptr<rayt::Camera> camera(std::make_unique<rayt::Camera>(x, y, z));

    for (int j = 0; j < ny; ++j) {
        std::cerr << "Rendering (y = " << j << ") " << (100.0 * j / (ny - 1)) << "%" << std::endl;
        for (int i = 0; i < nx; ++i) {
            float u = float(i) / float(nx);
            float v = float(j) / float(ny);
            rayt::Ray r = camera->getRay(u, v);
            vec3 c = color(r);
            image->write(i, j, c.getX(), c.getY(), c.getZ());
        }
    }
    stbi_write_bmp("render.bmp", nx, ny, sizeof(rayt::Image::rgb), image->pixels());

    return 0;
}