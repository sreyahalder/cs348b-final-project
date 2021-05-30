
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// lights/atmosphere.cpp*
#include "lights/atmosphere.h"
#include "imageio.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"
#include "math.h"

namespace pbrt {

Float hitAtmosphere(const Vector3f& center, Float radius, const Ray& r){
    Vector3f oc = Vector3f(r.o) - center;
    Float a = Dot(r.d, r.d);
    Float b = 2.0 * Dot(oc, r.d);
    Float c = Dot(oc, oc) - radius * radius;
    Float discriminant = b * b - 4 * a * c;
    if(discriminant < 0){
      return -1.0;
    } else {
      Float t0 = (-b + std::sqrt(discriminant)) / (2.0 * a);
      Float t1 = (-b - std::sqrt(discriminant)) / (2.0 * a);
      return t0 > 0 ? t0 : t1;
    }
}

// Atmosphere Method Definitions
Atmosphere::Atmosphere(const Spectrum &sunIntensity, const Vector3f &sunDir)
    : sunIntensity(sunIntensity),
      sunDir(sunDir) {
    betaRayleigh = Vector3f(3.8e-6f, 13.5e-6f, 33.1e-6f);
    betaMie = Vector3f(21e-6, 21e-6, 21e-6);
    betaAbsorption = Vector3f(2.04e-5, 4.97e-5, 1.95e-6);
    g = 0.76;
    Hr = 8500;
    Hm = 1200;
    absorptionHeightMax = 30000;
    absorptionFalloff = 3000;
    step_count = 16;
    light_step_count = 8;
    earthRadius = 6360e3;
    atmosphereRadius = 6420e3;
}

Vector3f Atmosphere::SampleLightRay(const Ray &ray) {
  Vector3f center(0., -earthRadius, 0.);
  Float hitT = hitAtmosphere(center, atmosphereRadius, ray);
  Point3f p = ray(hitT);
  Float stepDist = Distance(p, ray.o) / light_step_count;
  Float stepT = 0;

  Vector3f opticalDepthLight(0., 0., 0.);

  for (int j = 0; j < light_step_count; j++) {
    Point3f curPosition = ray(stepT + stepDist * 0.5f);
    Float h = curPosition.y;

    Float curOpticalDepthRayleigh = exp(-h / Hr) * stepDist;
    Float curOpticalDepthMie = exp(-h / Hm) * stepDist;
    Float curOpticalDepthOzone = 1.0 / cosh((absorptionHeightMax - h) / absorptionFalloff);
    curOpticalDepthOzone *= curOpticalDepthRayleigh * stepDist;

    opticalDepthLight += Vector3f(curOpticalDepthRayleigh, curOpticalDepthMie, curOpticalDepthOzone);

    stepT += stepDist;
  }
  return opticalDepthLight;
}

Spectrum Atmosphere::ComputeScattering(const Ray &ray, const SurfaceInteraction &isect, const bool renderAtmosphere) {
    Float mu = Dot(ray.d, sunDir);
    Float mumu = mu * mu;
    Float gg = g * g;
    Float phaseRayleigh = 3.0 / (16.0 * Pi) * (1.0 + mumu);
    Float phaseMie = 3.0 / (8.0 * Pi) * ((1.0 - gg) * (mumu + 1.0)) /
                     (std::pow(1.0 + gg - 2.0 * mu * g, 1.5) * (2.0 + gg));

    Point3f p = isect.p;
    if (renderAtmosphere) {
      Vector3f center(0., -earthRadius, 0.);
      Float hitT = hitAtmosphere(center, atmosphereRadius, ray);
      p = ray(hitT);
    }

    // calculate step size based on distance to surface intersection?
    Float stepDist = Distance(p, ray.o) / step_count;

    Vector3f accumulatedRayleigh(0., 0., 0.);
    Vector3f accumulatedMie(0., 0., 0.);
    Vector3f opticalDepth(0., 0., 0.);

    // start at t = 0
    Float stepT = 0;

    // loop over step count:
    for (int i = 0; i < step_count; i++) {
        // get point on ray at step t
        Point3f curPosition = ray(stepT + stepDist * 0.5f);

        // get height (altitude?) of point
        Float h = curPosition.y;

        // calculate current optical depth Rayleigh and Mie
        Float curOpticalDepthRayleigh = exp(-h / Hr) * stepDist;
        Float curOpticalDepthMie = exp(-h / Hm) * stepDist;

        // calculate current Optical Depth Ozone
        Float curOpticalDepthOzone = 1.0 / cosh((absorptionHeightMax - h) / absorptionFalloff);
        curOpticalDepthOzone *= curOpticalDepthRayleigh * stepDist;

        // add optical depth R M O to overall optical depth
        opticalDepth += Vector3f(curOpticalDepthRayleigh, curOpticalDepthMie, curOpticalDepthOzone);

        // Sample light ray at current ray step sample -> optical depth light
        Ray rayToSun(curPosition, sunDir);
        Vector3f opticalDepthLight = SampleLightRay(rayToSun);
        // Vector3f opticalDepthLight(1., 1., 1.);

        // calculate attenuation from betaR,M,A , optical depth and optical depth light
        Vector3f r = betaRayleigh * (opticalDepth.x + opticalDepthLight.x) + 
                     betaMie * 1.1f * (opticalDepth.y + opticalDepthLight.y);
        Vector3f attn(exp(-r.x), exp(-r.y), exp(-r.z));

        // calculate accumulated Rayleigh and Mie
        accumulatedRayleigh += curOpticalDepthRayleigh * attn;
        accumulatedMie += curOpticalDepthMie * attn;

        // increment primary step position
        stepT += stepDist;
    }

    // calculate scattering color based on sun intensity * phase functions * beta * accumulated rayleigh
    Spectrum scatteringColor = sunIntensity;
    scatteringColor[0] *= phaseRayleigh * betaRayleigh.x * accumulatedRayleigh.x + phaseMie * betaMie.x * accumulatedMie.x;
    scatteringColor[1] *= phaseRayleigh * betaRayleigh.y * accumulatedRayleigh.y + phaseMie * betaMie.y * accumulatedMie.y;
    scatteringColor[2] *= phaseRayleigh * betaRayleigh.z * accumulatedRayleigh.z + phaseMie * betaMie.z * accumulatedMie.z;

    if (isnan(scatteringColor[0]) || isnan(scatteringColor[1]) || isnan(scatteringColor[2])) {
      // return Spectrum(0.f);
      Float space[3] = {.036f, .056f, .102f};
      return Spectrum::FromRGB(space);
    }
    if (!renderAtmosphere) scatteringColor *= sunDir.y * 100;
    return scatteringColor;
}

std::shared_ptr<Atmosphere> CreateAtmosphere() {
  Spectrum sunIntensity(40.f);
  // Vector3f sunDir(0., 1., 0.);
  Vector3f sunDir(0., 0.01, -1.); //<-- Sunset
  // Spectrum sunIntensity(.5f);
  // Vector3f sunDir(0., 0.3, -1.);  //<-- Moonlight
  return std::make_shared<Atmosphere>(sunIntensity, Normalize(sunDir));
}

}  // namespace pbrt
