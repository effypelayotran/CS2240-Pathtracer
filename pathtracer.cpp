#include "pathtracer.h"

#include <iostream>
#include <cmath>
#include <random>
#include <tuple>
#include <scene/scene.h>

#include <Eigen/Dense>

#include <util/CS123Common.h>

using namespace Eigen;

/*
 * Random Float Generator for Monte Carlo Estiamtors
 */
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

/*
 * Diffuse BRDF Equation
 */
Vector3f diffuseBRDF(Vector3f matDiff) {
    return matDiff / M_PI;
}

/*
 * Glossy Specular BRDF Equation
 */
Vector3f glossyBRDF(Vector3f matSpec, float n, Vector3f wi, Vector3f wo, Vector3f normal) {
    Vector3f reflected_wi = wi - 2 *(wi.dot(normal))*normal;
    reflected_wi.normalize();
    return matSpec * ((n + 2.0) / (2.0 * M_PI)) * pow(reflected_wi.dot(wo), n);
}

/*
 * Reflection BRDF Equation (assumed w points into surface)
 */
Vector3f mirrorBRDF(Vector3f w, Vector3f normal) {
    return w - 2 *(w.dot(normal))*normal;
}


/*
 * Uniform Hemisphere Sampling
 */
Vector3f sampleNextDirection(const Vector3f &normal) {
    float e1 = dis(gen);
    float e2 = dis(gen);
    float phi = 2 * M_PI * e1;
    float theta = acos(e2);

    // Converting randomly selected Spherical Coordinates (phi & theta) into Cartesian Coordinates (xyz)
    Vector3f sample(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
        );

    // Creating an orthonormal basis about the normal
    Vector3f up;
    if (std::abs(normal.x()) < std::abs(normal.y()) && std::abs(normal.x()) < std::abs(normal.z())) {
        up = Vector3f(1, 0, 0);  // Choose X if X is the smallest component
    } else if (std::abs(normal.y()) < std::abs(normal.z())) {
        up = Vector3f(0, 1, 0);  // Choose Y if Y is smaller than Z
    } else {
        up = Vector3f(0, 0, 1);  // Otherwise, choose Z
    }
    Vector3f tangent = up.cross(normal).normalized();
    Vector3f bitangent = normal.cross(tangent);

    // Align randomly selected (xyz) vector about the normal by projecting it onto the orthonormal basis
    Vector3f alignedSample = (sample.x() * tangent + sample.y() * bitangent + sample.z() * normal).normalized();
    return alignedSample;
}

/*
 * Importance Sampling for Diffuse BRDF where pdf is proportional to cos(theta)
 */
std::tuple<float, Vector3f> cosineWeightedSample(const Vector3f& normal) {
    float e1 = dis(gen);
    float e2 = dis(gen);
    float phi = 2 * M_PI * e1;
    float theta = acos(e2);

    // Converting randomly selected Spherical Coordinates (phi & theta) into Cartesian Coordinates (xyz)
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    Vector3f sample(x, y, z);

    // Creating an orthonormal basis about the normal
    Vector3f up;
    if (std::abs(normal.x()) < std::abs(normal.y()) && std::abs(normal.x()) < std::abs(normal.z())) {
        up = Vector3f(1, 0, 0);
    } else if (std::abs(normal.y()) < std::abs(normal.z())) {
        up = Vector3f(0, 1, 0);
    } else {
        up = Vector3f(0, 0, 1);
    }
    Vector3f tangent = up.cross(normal).normalized();
    Vector3f bitangent = normal.cross(tangent);

    // Align randomly selected (xyz) vector about the normal by projecting it onto the orthonormal basis
    Vector3f alignedSample = (sample.x() * tangent + sample.y() * bitangent + sample.z() * normal).normalized();
    float pdf = cos(theta) / M_PI;
    return std::make_tuple(pdf, alignedSample);
}


/*
 * Importance Sampling for Glossy Specular BRDF where pdf is proportional to cos^n(theta)
 */
std::tuple<float, Vector3f> cosinePowerNWeightedSample(const Vector3f& normal, const Vector3f& wo, float n) {
    Eigen::Vector3f wo_refl = wo - 2 * wo.dot(normal) * normal;
    float e1 = dis(gen);
    float e2 = dis(gen);
    float phi = 2 * M_PI * e1;
    float theta = acos(pow(e2, 1 / (n + 1)));

    // Converting randomly selected Spherical Coordinates (phi & theta) into Cartesian Coordinates (xyz)
    float x = sin(theta) * cos(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(theta);
    Vector3f sample(x, y, z);

    // Creating an orthonormal basis about the normal
    Vector3f up;
    if (std::abs(normal.x()) < std::abs(normal.y()) && std::abs(normal.x()) < std::abs(normal.z())) {
        up = Vector3f(1, 0, 0);
    } else if (std::abs(normal.y()) < std::abs(normal.z())) {
        up = Vector3f(0, 1, 0);
    } else {
        up = Vector3f(0, 0, 1);
    }
    Vector3f tangent = up.cross(wo_refl).normalized();
    Vector3f bitangent = wo_refl.cross(tangent);

    // Align randomly selected (xyz) vector about the normal by projecting it onto the orthonormal basis
    Vector3f alignedSample = (sample.x() * tangent + sample.y() * bitangent + wo_refl * sample.z()).normalized();
    float pdf = (n + 1) * std::pow(cos(theta), n) / (2 * M_PI);
    return std::make_tuple(pdf, alignedSample);
}


/*
 * Refraction BRDF Probability
 *      Schlick's Approximation for Fresnel Equation that
 *      states as the Angle from normal (x-axis of graph) gets bigger, the Reflectance % gets bigger (y-axis of graph)
 *      at a rate of a positive, 5th degree polynomial (not exponential actually!)
 *      Note that this property is only for Insulating Materials.
 * */
float schlickApprox(float ni, float nt, Vector3f w_Incident, Vector3f normal) {
    float cosThetaIncident = normal.dot(-w_Incident);
    float r0 = pow((ni-nt) / (ni+nt), 2);
    float rTheta = r0 + (1-r0) * pow(1-cosThetaIncident, 5);
    return rTheta;
}


/*
 * Radiance function for Direct Lighting Integral using Monte Carlo Estimator
 * */
std::tuple<Vector3f, Vector3f> directL(Vector3f p, Vector3f n_p, const std::vector<Triangle*>& emissiveTris, const Scene& scene, int numPointsToSample) {
    // returns all triangles in the scene whose material has non-zero emission
    // const std::vector<Triangle*>& getEmissives() const { return m_emissives; };
    Vector3f totalLight(0, 0, 0);
    Vector3f totalOmega(0, 0, 0);

    // Loop through all M Triangles in Scene that is a light source
    for (Triangle* triangleA : emissiveTris) {
        Vector3<Vector3f> vertices = triangleA->getVertices();
        std::vector<Vector3f> points;

        // Randomly Sample N Points in Triangle using Barycentric Coordinates
        points.reserve(numPointsToSample);
        for (int i = 0; i < numPointsToSample; ++i) {
            float alpha = dis(gen);
            float beta = dis(gen);

            // alpha + beta + gamma = 1
            if (alpha + beta > 1.0f) {
                alpha = 1.0f - alpha;
                beta = 1.0f - beta;
            }
            // gamma is implicitly defined as 1 - alpha - beta = gamma
            float gamma = 1.0f - alpha - beta;

            Eigen::Vector3f point = alpha * vertices[0] + beta * vertices[1] + gamma * vertices[2];
            points.push_back(point);
        }

        // Compute Triangle Area as 1/2 || a x b ||  where are a & b any two sides of triangle
        float area;
        Vector3f a = vertices[1] - vertices[0];
        Vector3f b = vertices[2] - vertices[0];
        Vector3f axb = a.cross(b);
        float parallelogram = axb.norm();
        area = (0.5f) * parallelogram;


        // Loop through the N points in the Triangle and estimate Radiance from Triangle
        const tinyobj::material_t& mat = triangleA->getMaterial();
        Vector3f emittedLight = Vector3f(mat.emission[0], mat.emission[1], mat.emission[2]);
        Vector3f accLight = Vector3f(0, 0, 0);
        Vector3f accOmega = Vector3f(0, 0, 0);

        for(Vector3f p_prime:points) {
            // Vector3f w_prime = (p - (p_prime)).normalized();
            Vector3f w = (p_prime - p).normalized();
            Vector3f w_prime = -w;

            if (w_prime.y() >= 0) {
                continue;
            }

            IntersectionInfo i;
            Ray ray(p, w);
            if (scene.getIntersection(ray, &i)) {
                const Triangle *t = static_cast<const Triangle *>(i.data);
                const tinyobj::material_t& mat = t->getMaterial();

                // V(p, p') binary visibility function
                if (t->getIndex() == triangleA->getIndex()){
                    float cosTheta = w.dot(n_p);
                    float cosThetaPrime = w_prime.dot(triangleA->getNormal(i));
                    // cos(0)cos(0') / |p-p'|^2
                    float term = (cosTheta * cosThetaPrime) / (p - p_prime).squaredNorm();
                    // L(p', w') emittedLight from Area Light Source  A'
                    Vector3f L = emittedLight * term;
                    accLight += L;
                    accOmega += w;
                }
            }
        }

        Vector3f outputLight = accLight * (area / numPointsToSample);
        Vector3f outputOmega = (accOmega / numPointsToSample);
        totalLight += outputLight;
        totalOmega += outputOmega;
    }

    Vector3f averageLight = totalLight / emissiveTris.size();
    Vector3f averageOmega = totalOmega / emissiveTris.size();

    return std::make_tuple(averageLight, averageOmega);
}



PathTracer::PathTracer(int width, int height)
    : m_width(width), m_height(height)
{
}

void PathTracer::traceScene(QRgb *imageData, const Scene& scene)
{
    std::cout << settings.isLowDiscrepancySampling << std::endl;
    std::vector<Vector3f> intensityValues(m_width * m_height);
    Matrix4f invViewMat = (scene.getCamera().getScaleMatrix() * scene.getCamera().getViewMatrix()).inverse();
    for(int y = 0; y < m_height; ++y) {
        std::cout << "y= " << y << std::endl;
        #pragma omp parallel for
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            // intensityValues[offset] = tracePixel(x, y, scene, invViewMat);

            if (settings.isStratifiedSampling) {
                // Extra Feature: Stratified Sampling
                intensityValues[offset] = traceStratified(x, y, scene, invViewMat);
            } else if (settings.isLowDiscrepancySampling) {
                // Extra Feature: Low Discrepancy Sampling
                intensityValues[offset] = traceLowDiscrep(x, y, scene, invViewMat);
            } else {
                intensityValues[offset] = tracePixel(x, y, scene, invViewMat);
            }
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

        output += traceRay(r, scene, true);
    }
    return output / float(numSamples);
}


/*
 * Stratified Pixel Sampling
 * */
Vector3f PathTracer::traceStratified(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix)
{
    Vector3f output(0, 0, 0);
    int numSamples = settings.samplesPerPixel;

    // Create a grid with n squares total (sqrt(n) x sqrt(n) = n total squares) and sample 1 pizel per square
    int sqrtNumSamples = int(sqrt(numSamples));
    // numSamples taken in this pixel gets rounded down to a perfect square
    int recalculatedNumSamples = pow(sqrtNumSamples, 2);

    for (int i = 0; i < sqrtNumSamples; ++i) {
        for (int j = 0; j < sqrtNumSamples; ++j) {
            float new_x = (x + ((i + dis(gen)) / sqrtNumSamples)) / m_width;
            float new_y = (y + ((j + dis(gen)) / sqrtNumSamples)) / m_height;
            Vector3f p(0, 0, 0);
            Vector3f d(2.f * new_x - 1.0f, 1.0f - 2.f * new_y , -1);
            d.normalize();

            Ray r(p, d);
            r = r.transform(invViewMatrix);
            output += traceRay(r, scene, true);
        }
    }

    return output / recalculatedNumSamples;
}

/*
 * Van Der Corput Sequence for less clumpy random number generator
 * */
float vanDerCorput(int map_this_number, int base) {
    float result = 0.0f;
    float f = 1.0f / base;
    while (map_this_number > 0) {
        result += (map_this_number % base) * f;
        map_this_number /= base;
        f /= base;
    }
    return result;
}

/*
 * Low Discrepancy Pixel Sampling using Quasi-Monte Carlo
 * */
Vector3f PathTracer::traceLowDiscrep(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix) {
    Vector3f output = Vector3f(0,0,0);
    int numSamples = settings.samplesPerPixel;

    for (int i = 1; i < (numSamples + 1); i++) {
        float new_x = ((float)x + vanDerCorput(i, 2)) / m_width;
        float new_y = ((float)y + vanDerCorput(i, 3)) / m_height;

        Vector3f p(0, 0, 0);
        Vector3f d((2.f * new_x) - 1, 1 - (2.f * new_y), -1);
        d.normalize();

        Ray r(p, d);
        r = r.transform(invViewMatrix);
        output += traceRay(r, scene, true);
    }
    return output / numSamples;
}

Vector3f PathTracer::traceRay(const Ray& r, const Scene& scene, bool countEmitted)
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

    IntersectionInfo i;
    Ray ray(r); // Ray.d points away from the source p and towards the surface
    if (scene.getIntersection(ray, &i)) {
        const Triangle *t = static_cast<const Triangle *>(i.data);
        const tinyobj::material_t& mat = t->getMaterial();
        Vector3f emittedLight(mat.emission[0], mat.emission[1], mat.emission[2]);
        Vector3f normal = t->getNormal(i);
        Vector3f diffuseColor(mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);
        Vector3f specularColor(mat.specular[0], mat.specular[1], mat.specular[2]);
        Vector3f ambientColor(mat.ambient[0], mat.ambient[1], mat.ambient[2]);
        // Vector3f L = emittedLight;
        // From Slack: Treat kd of specular & shiny materials as 1, 1, 1 for milestone
        // if (diffuseColor.x() < 0.1f && diffuseColor.y() < 0.1f && diffuseColor.z() < 0.1f) {
        //     diffuseColor = Vector3f(1, 1, 1);
        // }
        // Intiial Test: Just Print out the Diffuse
        // Vector3f diffuseColor(mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);
        // return diffuseColor;

        // Direct Lighting (Lo,d)
        Vector3f L = Vector3f(0,0,0);
        auto directLighting = directL(i.hit, t->getNormal(i), scene.getEmissives(), scene, settings.numDirectLightingSamples);
        Vector3f directRadiance = std::get<0>(directLighting);
        Vector3f directDirection = std::get<1>(directLighting);
        if (mat.shininess < 500 && mat.ior <= 2) {
            Vector3f brdf_direct;
            if (specularColor.x() < 0.5) {
                brdf_direct = diffuseBRDF(diffuseColor);
            } else {
                brdf_direct = glossyBRDF(specularColor, mat.shininess, -directDirection, -r.d, t->getNormal(i));
            }
            L = directRadiance.cwiseProduct(brdf_direct);
        }

        // Indirect Lighting (Lo,i)
        if (!settings.directLightingOnly) {
            float pdf_rr = settings.pathContinuationProb;
            if (dis(gen) < pdf_rr) {
                // Diffuse & Glossy Specular BRDF
                if (mat.shininess < 500 && mat.ior <= 2) {
                    Vector3f incomingRadiance;
                    Vector3f wi;
                    float pdf;
                    Vector3f brdf;

                    if (specularColor.x() >= 0.5) {

                        if (!settings.isImportance) {
                            wi = sampleNextDirection(normal);
                            pdf = 1 / (2 * M_PI);
                        } else {
                            // Extra Feature: Importance Sampling for Glossy BRDF
                            auto sampled = cosinePowerNWeightedSample(normal, r.d, mat.shininess);
                            pdf = std::get<0>(sampled);
                            wi = std::get<1>(sampled);
                        }

                        Ray incomingRay(i.hit, wi);
                        incomingRadiance = traceRay(incomingRay, scene, false);
                        // For this BRDF calculation, wi should point into surface & wo should point away from surface hence the - signs here:
                        brdf = glossyBRDF(specularColor, mat.shininess, -incomingRay.d, -r.d, normal);
                    } else {
                        if (!settings.isImportance) {
                            wi = sampleNextDirection(normal);
                            pdf = 1 / (2 * M_PI);
                        } else {
                            // Extra Feature: Importance Sampling for Diffuse BRDF
                            auto sampled = cosineWeightedSample(normal);
                            pdf = std::get<0>(sampled);
                            wi = std::get<1>(sampled);
                        }

                        Ray incomingRay(i.hit, wi);
                        incomingRadiance = traceRay(incomingRay, scene, false);
                        brdf = diffuseBRDF(diffuseColor);
                    }

                    L += incomingRadiance.cwiseProduct(brdf) * wi.dot(normal) / (pdf * pdf_rr);
                }

                // Mirror Reflection BRDF
                if (mat.shininess >= 500) {
                    //For this BRDF, don't sample a random wi, just pick wi as the reflection ray.
                    //Vector3f brdf = (r.d - 2 * (r.d.dot(normal) * normal)).normalized();
                    Vector3f brdf = mirrorBRDF(r.d, normal);
                    Ray incomingRay(i.hit, brdf);
                    Vector3f incomingRadiance = traceRay(incomingRay, scene, true);
                    L += incomingRadiance / pdf_rr;
                }

                // Refraction with Fresnel Reflection BSDF
                if (mat.ior > 2.0) {
                    Vector3f wi;
                    float cosThetaIncident = r.d.dot(normal);

                    if (cosThetaIncident < 0) {
                        // 7a: If cosThetaIncident is negative, then wo entering the material from air (eta = 1.0)
                        float rTheta = schlickApprox(1.0, mat.ior, r.d, normal);
                        if (dis(gen) < rTheta) {
                            wi = r.d - 2 * r.d.dot(normal) * normal;
                        }

                        else {
                            float etaRatio = 1.0 / mat.ior;
                            float cosThetaIncident = normal.dot(-r.d);
                            float criticalcosTheta = (1 - etaRatio * etaRatio * (1 - cosThetaIncident * cosThetaIncident));
                            // From Snell, we solved that for cosThetaT = sqrt((1 - etaRatio * etaRatio * (1 - cosThetaIn * cosThetaIn)))
                            // But this term under the sqrt must be non-negative for a real solution to exist.
                            // If it's negative, the square root would involve an imaginary number,
                            // which physically means refraction cannot occur.

                            // Thus, TIR occurs when the term under the sqrt for
                            // cosTheta = sqrt((1 - etaRatio * etaRatio * (1 - cosThetaIn * cosThetaIn)))
                            // is negative, because the refracted angle thetaT could not exist in real number.

                            if (criticalcosTheta < 0) {
                                // 7b: Optical Manhole: Total Internal Reflection past the Critical Angle
                                wi = (r.d - 2 * r.d.dot(normal) * normal).normalized();
                            }
                            else {
                                // Refraction: GLSL refract function used safely within Critical Angle
                                float cosThetaOut = sqrt(1 - etaRatio * etaRatio * (1 - cosThetaIncident * cosThetaIncident));
                                wi = (etaRatio * r.d + (etaRatio*cosThetaIncident - cosThetaOut) * normal).normalized();
                            }
                        }
                    }

                    else {
                        // 7a: If cosThetaIncidnet is postive, then wo is exiting the material from air
                        float rTheta = schlickApprox(mat.ior, 1.0, r.d, -normal);
                        if (dis(gen) < rTheta) {
                            wi = r.d - 2 * r.d.dot(-normal) * -normal;
                        }

                        else {
                            float etaRatio = mat.ior / 1.0;
                            float cosThetaIn = (-normal).dot(-r.d);
                            float criticalcosTheta = (1 - etaRatio * etaRatio * (1 - cosThetaIn * cosThetaIn));
                            // From Snell, we solved that for cosThetaT = sqrt((1 - etaRatio * etaRatio * (1 - cosThetaIn * cosThetaIn)))
                            // But this term under the sqrt must be non-negative for a real solution to exist.
                            // If it's negative, the square root would involve an imaginary number,
                            // which physically means refraction cannot occur.

                            // Thus, TIR occurs when the term under the sqrt for
                            // cosTheta = sqrt((1 - etaRatio * etaRatio * (1 - cosThetaIn * cosThetaIn)))
                            // is negative, because the refracted angle thetaT could not exist in real number.

                            if (criticalcosTheta < 0) {
                                // 7b: Optical Manhole: Total Internal Reflection past the Critical Angle
                                wi = (r.d - 2 * r.d.dot(-normal) * -normal).normalized();
                            }
                            else {
                                // Refraction: GLSL refract function used safely within Critical Angle
                                float cosThetaT = sqrt(1 - etaRatio * etaRatio * (1 - cosThetaIn * cosThetaIn));
                                wi  = (etaRatio * r.d + (etaRatio * cosThetaIn - cosThetaT) * -normal).normalized();
                            }
                        }
                    }

                    // Wi points away from the surface by convention
                    Ray incomingRay(i.hit, wi);
                    float attenuation = 1.0;
                    // std::cout << settings.isAttenuate << std::endl;

                    // Extra Feature: Attenuate Refracted Light
                    if (settings.isAttenuate && (wi.dot(t->getNormal(i)) < 0)) {
                        Ray ray(i.hit, wi);
                        IntersectionInfo i2;
                        if (scene.getIntersection(ray, &i2)) {
                            float dist = (i.hit - i2.hit).norm();
                            attenuation = 1.0 / (1.0 + 1.0 * dist + 1.0* pow(dist, 2));
                            // std::cout << "Hi:)" << std::endl;
                        }
                    }

                    Vector3f incomingRadiance = attenuation * traceRay(incomingRay, scene, true);
                    L += incomingRadiance / pdf_rr;
                }
            }
        }


        if (countEmitted) {
            L += Vector3f(mat.emission[0], mat.emission[1], mat.emission[2]);
        }

        return L;
    } else {
        // Background color
        return Vector3f(0.0f, 0.0f, 0.0f);
    }
}

/*
 * Extended Reinhard TM Function where whitest white gets mapped to 1.0
 * */
float extendedReinhard(float C, float C_white) {
    return (C * (1 + (C / pow(C_white, 2)))) / (1 + C);
}

float reinhard(float L) {
    return (L / (1.0f + L));
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
            // Reinhard is a mapping function where L_mapped =((L / 1 + L)
            // If L < 1, then L / (1 + L) ≈ L, so dark colors stay dark.
            // If L > 1, then L / (1 + L) ≈ 1, preventing bright spots from blowing out.

            // High dynamic range color from path tracing
            Vector3f hdrColor = intensityValues[offset];

            // Luminance scaling
            float luminance = 0.2126f * hdrColor.x() + 0.7152f * hdrColor.y() + 0.0722f * hdrColor.z();
            float C_white = 0.8f;
            float gamma = 2.2f;

            // Normalize RGB values using Extended Reinhard Luminance ratio
            // float mappedLuminance = extendedReinhard(luminance, C_white);
            float mappedLuminance = reinhard(luminance);
            Vector3f ldrColor = hdrColor * (mappedLuminance / luminance);

            // Gamma correction
            ldrColor.x() = std::pow(ldrColor.x(), 1.0f / gamma);
            ldrColor.y() = std::pow(ldrColor.y(), 1.0f / gamma);
            ldrColor.z() = std::pow(ldrColor.z(), 1.0f / gamma);

            // Clamp and scale to [0, 255] 8-bit range
            int r = std::min(255, std::max(0, static_cast<int>(ldrColor.x() * 255)));
            int g = std::min(255, std::max(0, static_cast<int>(ldrColor.y() * 255)));
            int b = std::min(255, std::max(0, static_cast<int>(ldrColor.z() * 255)));

            imageData[offset] = qRgb(r, g, b);
        }
    }
}


