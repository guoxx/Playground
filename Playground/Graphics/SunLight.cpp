/***************************************************************************
# Copyright (c) 2015, NVIDIA CORPORATION. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of NVIDIA CORPORATION nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***************************************************************************/
#include "Framework.h"
#include "SunLight.h"
#include "Utils/Gui.h"
#include "API/Device.h"
#include "API/ConstantBuffer.h"
#include "API/Buffer.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "Data/VertexAttrib.h"
#include "Graphics/Model/Model.h"
#include "Utils/Spectrum.h"
#include "HosekWilkie_SkylightModel/ArHosekSkyModel.h"


namespace Falcor
{
    /* Apparent radius of the sun as seen from the earth (in degrees).
    This is an approximation--the actual value is somewhere between
    0.526 and 0.545 depending on the time of year */
    #define SUN_APP_ANGULAR_DIAMETER 0.5358f

    glm::vec2 ConcentricSampleDisk(const glm::vec2 &u)
    {
        // Map uniform random numbers to $[-1,1]^2$
        glm::vec2 uOffset = 2.f * u - glm::vec2(1, 1);

        // Handle degeneracy at the origin
        if (uOffset.x == 0 && uOffset.y == 0) return glm::vec2(0, 0);

        // Apply concentric mapping to point
        Float theta, r;
        if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
            r = uOffset.x;
            theta = float(M_PI_4 * (uOffset.y / uOffset.x));
        }
        else {
            r = uOffset.y;
            theta = float(M_PI_2 - M_PI_4 * (uOffset.x / uOffset.y));
        }
        return r * glm::vec2(std::cos(theta), std::sin(theta));
    }

	SunLight::SunLight()
    {
        mTurbidity = 2.0f;
        mTheta = glm::radians(45.0f);
        mPhi = glm::radians(20.0f);
        mGroundAlbedo = glm::vec3(0.5, 0.5, 0.5);

        updateLightInfo();
    }

    SunLight::SharedPtr SunLight::create()
    {
        SunLight* pLight = new SunLight();
        return SharedPtr(pLight);
    }

    SunLight::~SunLight() = default;

    void SunLight::renderUI(Gui* pGui, const char* group)
    {
        if(!group || pGui->beginGroup(group))
        {
            bool needToUpdate = false;
            if (pGui->addFloatVar("Theta", mTheta, 0, float(M_PI)))
            {
                needToUpdate = true;
            }
            if (pGui->addFloatVar("Phi", mPhi, 0, float(M_PI * 2)))
            {
                needToUpdate = true;
            }
            if (pGui->addFloatVar("Turbidity", mTurbidity, 1, 10, 1))
            {
                needToUpdate = true;
            }
            if (pGui->addRgbColor("GroundAlbedo", mGroundAlbedo))
            {
                needToUpdate = true;
            }
            Light::renderUI(pGui);
            if (group)
            {
                pGui->endGroup();
            }

            if (needToUpdate)
            {
                updateLightInfo();
            }
        }
    }

    void SunLight::prepareGPUData()
    {
    }

    void SunLight::unloadGPUData()
    {
    }

    void SunLight::move(const glm::vec3& position, const glm::vec3& target, const glm::vec3& up)
    {
        logError("SunLight::move() is not used and thus not implemented for now.");
    }

    SampledSpectrum SunLight::computeSunRadiance() const
    {
        const SphericalCoordinates sphericalCoord = SphericalCoordinates::FromThetaAndPhi(mTheta, mPhi);
        const SampledSpectrum groundAlbedoSpectrum = SampledSpectrum::FromRGB(mGroundAlbedo);

        SampledSpectrum sunRadiance;

        const uint64 NumDiscSamples = 8;
        for (uint64 x = 0; x < NumDiscSamples; ++x)
        {
            for (uint64 y = 0; y < NumDiscSamples; ++y)
            {
                float u = (x + 0.5f) / NumDiscSamples;
                float v = (y + 0.5f) / NumDiscSamples;
                glm::vec2 discSamplePos = ConcentricSampleDisk(glm::vec2(u, v));

                float elevation = sphericalCoord.GetElevation();
                float theta = sphericalCoord.GetZenith() + discSamplePos.y * glm::radians(SUN_APP_ANGULAR_DIAMETER * 0.5f);
                float gamma = discSamplePos.x * glm::radians(SUN_APP_ANGULAR_DIAMETER * 0.5f);

                SampledSpectrum solarRadiance;
                for (int32 i = 0; i < NumSpectralSamples; ++i)
                {
                    ArHosekSkyModelState* skyState = arhosekskymodelstate_alloc_init(elevation, mTurbidity, groundAlbedoSpectrum[i]);
                    float wavelength = glm::lerp(float(SampledLambdaStart), float(SampledLambdaEnd), i / float(NumSpectralSamples));

                    solarRadiance[i] = float(arhosekskymodel_solar_radiance(skyState, theta, gamma, wavelength));

                    arhosekskymodelstate_free(skyState);
                    skyState = nullptr;
                }

                sunRadiance += solarRadiance * (1.0f / NumDiscSamples) * (1.0f / NumDiscSamples);
            }
        }
        return sunRadiance;
    }

    void SunLight::updateLightInfo()
    {
        SampledSpectrum sunRadianceSPD = computeSunRadiance();

        double sunRadiance = 0;
        double sunLuminance = 0;
        for (int i = 0; i < NumSpectralSamples; ++i)
        {
            sunRadiance += sunRadianceSPD[i] * SpectrumSamplesStep;
            sunLuminance += sunRadianceSPD[i] * 683 * SampledSpectrum::Y[i] * SpectrumSamplesStep;
        }
        double luminousEfficiency = sunLuminance / sunRadiance;

        RGBSpectrum radiance = sunRadianceSPD.ToRGBSpectrum();
        RGBSpectrum luminance = radiance * CIE_Y_integral * float(luminousEfficiency);

        float theta = glm::radians(SUN_APP_ANGULAR_DIAMETER * 0.5f);
        float solidAngle = 2 * float(M_PI) * (1 - std::cos(theta));
        RGBSpectrum illuminance = luminance * solidAngle;

        float maxIlluminance = std::max(illuminance[0], std::max(illuminance[1], illuminance[2]));
        RGBSpectrum normalizedIlluminance = illuminance / maxIlluminance;
        setColorFromUI(normalizedIlluminance.ToRGB());
        setIntensityFromUI(maxIlluminance);

        const SphericalCoordinates sphericalCoord = SphericalCoordinates::FromThetaAndPhi(mTheta, mPhi);
        setWorldDirection(-SphericalCoordinates::ToSphere(sphericalCoord));
    }
}
