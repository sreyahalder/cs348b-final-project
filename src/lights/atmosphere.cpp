
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

// Atmosphere Method Definitions
Atmosphere::Atmosphere(const Spectrum &sunIntensity, const Vector3f &sunDir)
    : sunIntensity(sunIntensity),
      sunDir(sunDir) {
    betaRayleigh = Vector3f(5.5e-6, 13.0e-6, 22.4e-6);
    betaMie = Vector3f(21e-6, 21e-6, 21e-6);
    betaAbsorption = Vector3f(2.04e-5, 4.97e-5, 1.95e-6);
    g = 0.76;
    rayleighScale = 8500;
    mieScale = 1200;
    absorptionHeightMax = 30000;
    absorptionFalloff = 3000;
}

Spectrum Atmosphere::ComputeScattering(const Ray &ray, const SurfaceInteraction &isect) {
    Float mu = Dot(ray.d, sunDir);
    Float mumu = mu * mu;
    Float gg = g * g;
    Float phaseRayleigh = 3.0 / (16.0 * Pi) * (1.0 + mumu);
    Float phaseMie = 3.0 / (8.0 * Pi) * ((1.0 - gg) * (mumu + 1.0)) / (std::pow(1.0 + gg - 2.0 * mu * g, 1.5) * (2.0 + gg));

    // calculate step size based on distance to surface intersection?
    Float STEP_COUNT = 5;
    Float stepDist = Distance(isect.p, ray.o) / STEP_COUNT;
    Float stepSize = stepDist;

    Vector3f accumulatedRayleigh(0., 0., 0.);
    Vector3f accumulatedMie(0., 0., 0.);
    Vector3f opticalDepth(0., 0., 0.);

    // start at t = 0
    Float stepT = stepSize;

    // loop over step count:
    for (int i = 0; i < STEP_COUNT; i++) {
        // get point on ray at step t
        Point3f curPosition = ray(stepT);

        // get height (altitude?) of point
        Float curHeight = curPosition.y;

        // calculate current optical depth Rayleigh and Mie
        Float curOpticalDepthRayleigh = exp(-curHeight / rayleighScale) * stepDist;
        Float curOpticalDepthMie = exp(-curHeight / mieScale) * stepDist;

        // calculate current Optical Depth Ozone
        Float curOpticalDepthOzone = 1.0 / cosh((absorptionHeightMax - curHeight) / absorptionFalloff);
        curOpticalDepthOzone *= curOpticalDepthRayleigh * stepDist;

        // add optical depth R M O to overall optical depth
        opticalDepth += Vector3f(curOpticalDepthRayleigh, curOpticalDepthMie, curOpticalDepthOzone);

        // Sample light ray at current ray step sample -> optical depth light
        // Vector3f opticalDepthLight = SampleLightRay(curPosition);
        Vector3f opticalDepthLight(1., 1., 1.); // <- abitrarily chosen

        // QUESTION: when we sample light rays, do we find the intersection between our step
        // along the original ray and the upper atmosphere? 

        // calculate attenuation from betaR,M,A , optical depth and optical depth light
        Vector3f r = betaRayleigh * (opticalDepth.x + opticalDepthLight.x) + 
                     betaMie * (opticalDepth.y + opticalDepthLight.y) +
                     betaAbsorption * (opticalDepth.z + opticalDepthLight.z);
        Vector3f attn(exp(-r.x), exp(-r.y), exp(-r.z));

        // calculate accumulated Rayleigh and Mie
        accumulatedRayleigh += curOpticalDepthRayleigh * attn;
        accumulatedMie += curOpticalDepthMie * attn;

        // increment primary step position
        stepT += stepSize;
    }

    // calculate scattering color based on sun intensity * phase functions * beta * accumulated rayleigh
    Vector3f R(phaseRayleigh * betaRayleigh.x * accumulatedRayleigh.x,
               phaseRayleigh * betaRayleigh.y * accumulatedRayleigh.y,
               phaseRayleigh * betaRayleigh.z * accumulatedRayleigh.z);
    Vector3f M(phaseMie * betaMie.x * accumulatedMie.x,
               phaseMie * betaMie.y * accumulatedMie.y,
               phaseMie * betaMie.z * accumulatedMie.z);
    Spectrum scatteringColor;
    scatteringColor[0] = sunIntensity[0] * (R + M).x;
    scatteringColor[1] = sunIntensity[1] * (R + M).y;
    scatteringColor[2] = sunIntensity[2] * (R + M).z;

    Vector3f sO = betaMie * opticalDepth.y + betaRayleigh * opticalDepth.x + betaAbsorption * opticalDepth.z;

    // calculate opacity of color based on same ^ + optical depth z
    Vector3f scatteringOpacity(exp(-sO.x), exp(-sO.y), exp(-sO.z));
    scatteringColor[0] *= scatteringOpacity.x; 
    scatteringColor[1] *= scatteringOpacity.y;
    scatteringColor[2] *= scatteringOpacity.z;

    return scatteringColor * 100;

    // Questions: how do we get the actual atmosphere background? Since it exists at infinity?
    // or do we just have a max depth for our ray and just calculate the color for that pixel?
    // How do we know what is the right distance for atmospheric scattering?
    // Can we use the same model as above to compute scattering color and opacity?
    // So since we are at max optical depth, would opacity just be 1? and essentially ignored?
}

std::shared_ptr<Atmosphere> CreateAtmosphere() {
    Spectrum sunIntensity(40.0f);
    Vector3f sunDir(0., 1., -1.);
    return std::make_shared<Atmosphere>(sunIntensity, Normalize(sunDir));
}

}  // namespace pbrt
