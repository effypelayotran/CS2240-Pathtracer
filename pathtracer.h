#ifndef PATHTRACER_H
#define PATHTRACER_H

#include <QImage>
#include "scene/scene.h"


struct Settings {
    int samplesPerPixel;
    bool directLightingOnly; // if true, ignore indirect lighting
    int numDirectLightingSamples; // number of shadow rays to trace from each intersection point
    float pathContinuationProb; // probability of spawning a new secondary ray == (1-pathTerminationProb)

    bool isStratifiedSampling = false;
    bool isLowDiscrepancySampling = false;
    bool isImportance = false;
    bool isAttenuate = false;
};

class PathTracer
{
public:
    PathTracer(int width, int height);

    void traceScene(QRgb *imageData, const Scene &scene);
    Settings settings;
    std::string pfmFilename;

private:
    int m_width, m_height;

    void toneMap(QRgb *imageData, std::vector<Eigen::Vector3f> &intensityValues);

    Vector3f tracePixel(int x, int y, const Scene &scene, const Eigen::Matrix4f &invViewMatrix);
    Vector3f traceLowDiscrep(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix);
    Vector3f traceStratified(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix);
    Vector3f traceRay(const Ray& r, const Scene &scene, bool countEmitted);
};

#endif // PATHTRACER_H
