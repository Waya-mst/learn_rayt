#include "rayt.h"

namespace rayt {

    class Shape;
    typedef std::shared_ptr<Shape> ShapePtr;

    class HitRec {
    public:
        float t;
        vec3 p;
        vec3 n;
    };


    class Shape {
    public:
        virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const = 0;
    };

    class ShapeList : public Shape {
    public:
        ShapeList() {}

        void add(const ShapePtr& shape) {
            m_list.push_back(shape);
        }

        virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
            HitRec temp_rec;
            bool hit_anything = false;
            float closest_so_far = t1;
            for (auto& p : m_list) {
                if (p->hit(r, t0, closest_so_far, temp_rec)) {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    hrec = temp_rec;
                }
            }

            return hit_anything;
        }

    private:
        std::vector<ShapePtr> m_list;
    };


    class Sphere : public Shape {
    public:
        Sphere() {}
        Sphere(const vec3& c, float r)
            :m_center(c)
            , m_radius(r) {}

        virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
            vec3 oc = r.origin() - m_center;
            float a = dot(r.direction(), r.direction());
            float b = 2.0f * dot(oc, r.direction());
            float c = dot(oc, oc) - pow2(m_radius);
            float D = b * b - 4 * a * c;

            if (D > 0) {
                float root = sqrtf(D);
                float temp = (-b - root) / (2.0f * a);
                if (temp < t1 && temp > t0) {
                    hrec.t = temp;
                    hrec.p = r.at(hrec.t);
                    hrec.n = (hrec.p - m_center) / m_radius;
                    return true;
                }
                temp = (-b + root) / (2.0f * a);
                if (temp < t1 && temp > t0) {
                    hrec.t = temp;
                    hrec.p = r.at(hrec.t);
                    hrec.n = (hrec.p - m_center) / m_radius;
                    return true;
                }
            }

            return false;
        }

    private:
        vec3 m_center;
        float m_radius;


    };

    class Scene {
    public:
        Scene(int width, int height)
            : m_image(std::make_unique<Image>(width, height))
            , m_backColor(0.2f) {}

        void build() {
            // camera

            vec3 w(-2.0f, -1.0f, -1.0f);
            vec3 u(4.0f, 0.0f, 0.0f);
            vec3 v(0.0f, 2.0f, 0.0f);
            m_camera = std::make_unique<Camera>(u, v, w);

            // Shapes

            ShapeList* world = new ShapeList();
            world->add(std::make_shared<Sphere>(
                vec3(0, 0, -1), 0.5f));
            world->add(std::make_shared<Sphere>(
                vec3(0, -100.5, -1), 100));
            m_world.reset(world);
        }


        vec3 color(const rayt::Ray& r, const Shape* world) const {
            HitRec hrec;
            vec3 c(0, 0, -1);

            if (world->hit(r, 0, FLT_MAX, hrec)) {
                return 0.5f * (hrec.n + vec3(1.0f));
            }
            return backgroundSky(r.direction());
        }

        vec3 background(const vec3& d) const {
            return m_backColor;
        }

        vec3 backgroundSky(const vec3& d) const {
            vec3 v = normalize(d);
            float t = 0.5f * (v.getY() + 1.0f);
            return lerp(t, vec3(1.0f), vec3(0.5f, 0.7f, 1.0f));
        }

        void render() {
            build();

            int nx = m_image->width();
            int ny = m_image->height();
#pragma omp parallel for schedule(dynamic, 1) num_threads(NUM_THREAD)
            for (int j = 0; j < ny; ++j) {
                std::cerr << "Rendering (y = " << j << ") " << (100.0 * j / (ny - 1)) << "%" << std::endl;
                for (int i = 0; i < nx; ++i) {

                    float u = float(i + drand48()) / float(nx);
                    float v = float(j + drand48()) / float(ny);
                    Ray r = m_camera->getRay(u, v);
                    vec3 c = color(r, m_world.get());
                    m_image->write(i, (ny - j - 1), c.getX(), c.getY(), c.getZ());
                }
            }

            stbi_write_bmp("render.bmp", nx, ny, sizeof(Image::rgb), m_image->pixels());
        }

    private:
        std::unique_ptr<Camera> m_camera;
        std::unique_ptr<Image> m_image;
        std::unique_ptr<Shape> m_world;
        vec3 m_backColor;
    };
}

int main()
{
    int nx = 200;
    int ny = 100;
    std::unique_ptr<rayt::Scene> scene(std::make_unique<rayt::Scene>(nx, ny));

    scene->render();

    return 0;
}