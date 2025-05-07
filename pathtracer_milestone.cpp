#include "pathtracer.h"

#include <iostream>
#include <random>
#include <scene/scene.h>

#include <Eigen/Dense>

#include <util/CS123Common.h>

using namespace Eigen;

// Init random float generator
std::random_device rd; // random seed from os
std::mt19937 gen(rd()); // pseudo-random generator
std::uniform_real_distribution<> dis(0.0, 1.0); // uniform distribution

Vector3f alignSampleToNormal(const Vector3f &sample, const Vector3f &normal) {
    Vector3f up(0, 1, 0);
    if (floatEpsEqual(normal.dot(up), 1.0f)) {
        up = Vector3f(1, 0, 0); // Handle edge case where normal is aligned with "up"/y
    }
    Vector3f tangent = up.cross(normal).normalized();
    Vector3f bitangent = normal.cross(tangent);

    return (sample.x() * tangent + sample.y() * bitangent + sample.z() * normal).normalized();
}

Vector3f sampleNextDirection(const Vector3f &normal) {
    float e1 = dis(gen);
    float e2 = dis(gen);
    float theta = acos(e1);
    float phi = 2.0f * M_PI * e2;

    // Convert random sample phi & theta into xyz cartesian coordinates
    Vector3f omega(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
        );

    // Align random sample xyz vector with the surface normal
    return alignSampleToNormal(omega, normal);
}


PathTracer::PathTracer(int width, int height)
    : m_width(width), m_height(height)
{
}

void PathTracer::traceScene(QRgb *imageData, const Scene& scene)
{
    std::vector<Vector3f> intensityValues(m_width * m_height);
    Matrix4f invViewMat = (scene.getCamera().getScaleMatrix() * scene.getCamera().getViewMatrix()).inverse();
    for(int y = 0; y < m_height; ++y) {
        std::cout << "y= " << y << std::endl;
        //#pragma omp parallel for
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            intensityValues[offset] = tracePixel(x, y, scene, invViewMat);
        }
    }

    toneMap(imageData, intensityValues);
}

Vector3f PathTracer::tracePixel(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix)
{
    Vector3f output(0, 0, 0);
    int numSamples = settings.samplesPerPixel;

    for (int s = 0; s < numSamples; ++s) {
        float px = (x + dis(gen)) / m_width;
        float py = (y + dis(gen)) / m_height;

        Vector3f p(0, 0, 0);
        Vector3f d(2.f * px - 1.f,
                   1.f - 2.f * py,
                   -1.f);
        d.normalize();

        Ray r(p, d);
        r = r.transform(invViewMatrix);

        output += traceRay(r, scene, 0);
    }

    return output / float(numSamples);
}


Vector3f PathTracer::traceRay(const Ray& r, const Scene& scene, int bounces)
{
    // IntersectionInfo i;
    // Ray ray(r);
    // if(scene.getIntersection(ray, &i)) {
    //     //** Example code for accessing materials provided by a .mtl file **
    //     //        const Triangle *t = static_cast<const Triangle *>(i.data);//Get the triangle in the mesh that was intersected
    //     //        const tinyobj::material_t& mat = t->getMaterial();//Get the material of the triangle from the mesh
    //     //        const tinyobj::real_t *d = mat.diffuse;//Diffuse color as array of floats
    //     //        const std::string diffuseTex = mat.diffuse_texname;//Diffuse texture name
    //     return Vector3f(1, 1, 1);
    // } else {
    //     return Vector3f(0, 0, 0);
    // }
    // bounces +=1;
    // std::cout << "bounce= " << bounces << std::endl;

    IntersectionInfo i;
    Ray ray(r);
    if (scene.getIntersection(ray, &i)) {
        // intiial test: return just the diffuse
        // const Triangle *t = static_cast<const Triangle *>(i.data);
        // const tinyobj::material_t& mat = t->getMaterial();
        // Vector3f diffuseColor(mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);
        // return diffuseColor;

        const Triangle *t = static_cast<const Triangle *>(i.data);
        const tinyobj::material_t& mat = t->getMaterial();
        Vector3f emittedLight(mat.emission[0], mat.emission[1], mat.emission[2]);
        Vector3f normal = t->getNormal(i);
        Vector3f diffuseColor(mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);
        Vector3f specularColor(mat.specular[0], mat.specular[1], mat.specular[2]);
        Vector3f ambientColor(mat.ambient[0], mat.ambient[1], mat.ambient[2]);


        Vector3f L = emittedLight;
        // From Slack: Treat kd of specular & shiny materials as 1, 1, 1 for milestone
        if (diffuseColor.x() < 0.1f && diffuseColor.y() < 0.1f && diffuseColor.z() < 0.1f) {
            diffuseColor = Vector3f(1, 1, 1);
        }

        // Russian Roulette termination
        float pdf_rr = 0.98;
        if (dis(gen) > pdf_rr) {
            return L;

        }
            // Uniform sampling
            Vector3f wi = sampleNextDirection(normal);
            float pdf = 1 / (2 * M_PI);
            // float pdf = wi.dot(normal) / M_PI;

            Ray incomingRay(i.hit + normal * 1e-5f, wi);
            Vector3f incomingRadiance = traceRay(incomingRay, scene, bounces);
            Vector3f brdf = diffuseColor / M_PI;

            // cwiseProduct = element wise multiplication
            L += incomingRadiance.cwiseProduct(brdf) * wi.dot(normal) / (pdf * pdf_rr);


            return L;

    } else {
        // Background color
        return Vector3f(0.0f, 0.0f, 0.0f);
    }
}

void PathTracer::toneMap(QRgb *imageData, std::vector<Vector3f> &intensityValues) {
    // for(int y = 0; y < m_height; ++y) {
    //     for(int x = 0; x < m_width; ++x) {
    //         int offset = x + (y * m_width);
    //         imageData[offset] = intensityValues[offset].norm() > 0 ? qRgb(255, 255, 255) : qRgb(40, 40, 40);
    //     }
    // }
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            // hdr = high dynamic range
            // ldr = low dynamic range
            // L)mapped =((L / 1 + L)
            // If L ≪ 1, then L / (1 + L) ≈ L, so dark colors stay dark.
            // If L ≫ 1, then L / (1 + L) ≈ 1, preventing bright spots from blowing out.
            Vector3f hdrColor = intensityValues[offset];
            Vector3f one_plus_L = Vector3f(1.0f + hdrColor.x(), 1.0f + hdrColor.y(), 1.0f + hdrColor.z());
            Vector3f ldrColor = hdrColor.cwiseQuotient(one_plus_L);

            // scale: L / (1 + L) is between [0, 1] multiply by 255 to get it between [0, 255]
            int r = static_cast<int>(ldrColor.x() * 255);
            int g = static_cast<int>(ldrColor.y() * 255);
            int b = static_cast<int>(ldrColor.x() * 255);

            // clamp: if it greater than 255 set to 255.
            r = std::min(255, r);
            g = std::min(255, g);
            b = std::min(255, b);

            imageData[offset] = qRgb(
                r,
                g,
                b
                );
        }
    }

}


