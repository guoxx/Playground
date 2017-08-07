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
#include "Utils/ThreadPool.h"


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

    double LuminousEfficiency(SampledSpectrum radiance)
    {
        double sunRadiance = 0;
        double sunLuminance = 0;
        for (int i = 0; i < NumSpectralSamples; ++i)
        {
            sunRadiance += radiance[i] * SpectrumSamplesStep;
            sunLuminance += radiance[i] * 683 * SampledSpectrum::Y[i] * SpectrumSamplesStep;
        }
        double luminousEfficiency = sunLuminance / sunRadiance;
        return luminousEfficiency;
    }

	SunLight::SunLight()
    {
        mData.mTurbidity = 2.0f;
        mData.mTheta = glm::radians(45.0f);
        mData.mPhi = glm::radians(20.0f);
        mData.mGroundAlbedo = glm::vec3(0.5, 0.5, 0.5);

        updateAsnyc(mData);
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
            if (pGui->addFloatVar("Theta", mData.mTheta, 0, float(M_PI_2) - 0.01f))
            {
                needToUpdate = true;
            }
            if (pGui->addFloatVar("Phi", mData.mPhi, 0, float(M_PI * 2)))
            {
                needToUpdate = true;
            }
            if (pGui->addFloatVar("Turbidity", mData.mTurbidity, 1.7f, 10.0f, 0.1f))
            {
                needToUpdate = true;
            }
            if (pGui->addRgbColor("GroundAlbedo", mData.mGroundAlbedo))
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
                updateAsnyc(mData);
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

    Texture::SharedPtr SunLight::GetSkyEnvMap() const
    {
        return mEnvMap;
    }

    SampledSpectrum SunLight::computeSunRadiance(float sunTheta, float sunPhi, float turbidity, glm::vec3 groundAlbedo)
    {
        const SphericalCoordinates sphericalCoord = SphericalCoordinates::FromThetaAndPhi(sunTheta, sunPhi);
        const SampledSpectrum groundAlbedoSpectrum = SampledSpectrum::FromRGB(groundAlbedo);

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
                    ArHosekSkyModelState* skyState = arhosekskymodelstate_alloc_init(elevation, turbidity, groundAlbedoSpectrum[i]);
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

    void SunLight::updateLightInfo(const InternalData& data)
    {
        SampledSpectrum sunRadianceSPD = computeSunRadiance(data.mTheta, data.mPhi, data.mTurbidity, data.mGroundAlbedo);

        double luminousEfficiency = LuminousEfficiency(sunRadianceSPD);

        RGBSpectrum radiance = sunRadianceSPD.ToRGBSpectrum();
        RGBSpectrum luminance = radiance * CIE_Y_integral * float(luminousEfficiency);

        float theta = glm::radians(SUN_APP_ANGULAR_DIAMETER * 0.5f);
        float solidAngle = 2 * float(M_PI) * (1 - std::cos(theta));
        RGBSpectrum illuminance = luminance * solidAngle;

        float maxIlluminance = std::max(illuminance[0], std::max(illuminance[1], illuminance[2]));
        RGBSpectrum normalizedIlluminance = illuminance / maxIlluminance;
        setColorFromUI(normalizedIlluminance.ToRGB());
        setIntensityFromUI(maxIlluminance);

        const SphericalCoordinates sphericalCoord = SphericalCoordinates::FromThetaAndPhi(data.mTheta, data.mPhi);
        setWorldDirection(-SphericalCoordinates::ToSphere(sphericalCoord));
    }

    void SunLight::updateEnvironmentMap(const InternalData& data)
    {
        const SphericalCoordinates sunCoord = SphericalCoordinates::FromThetaAndPhi(data.mTheta, data.mPhi);
        const vec3 sunDir = SphericalCoordinates::ToSphere(sunCoord);

        int32_t resolution = 1024;

        int32_t nTheta = resolution;
        int32_t nPhi = resolution * 2;

        float* pixelsData = new float[sizeof(float) * 4 * nPhi * nTheta];

        auto storePixel = [&](int32_t x, int32_t y, float r, float g, float b)
        {
            int32_t idx = (x + y * nPhi) * 4;
            pixelsData[idx] = r;
            pixelsData[idx + 1] = g;
            pixelsData[idx + 2] = b;
            pixelsData[idx + 3] = 1.0;
        };

        ArHosekSkyModelState* skyStates[3];

        for (int i = 0; i < _countof(skyStates); ++i)
        {
            skyStates[i] = arhosek_rgb_skymodelstate_alloc_init(data.mTurbidity, data.mGroundAlbedo[i], sunCoord.GetElevation());
        }

        for (int32_t t = 0; t < nTheta; ++t)
        {
            float theta = (t + 0.5f) / nTheta * float(M_PI);

            if (theta > M_PI_2)
            {
                for (int32_t p = 0; p < nPhi; ++p)
                {
                    storePixel(p, t, 0, 0, 0);
                }
                continue;
            }

            for (int32_t p = 0; p < nPhi; ++p)
            {
                float phi = (p + 0.5f) / nPhi * float(M_PI) * 2;

                SphericalCoordinates sphericalCoord = SphericalCoordinates::FromThetaAndPhi(theta, phi);
                glm::vec3 dir = SphericalCoordinates::ToSphere(sphericalCoord);

                float gamma = std::acos(glm::clamp(glm::dot(dir, sunDir), -1.0f, 1.0f));

                RGBSpectrum radiance;
                for (int i = 0; i < _countof(skyStates); ++i)
                {
                    radiance[i] = float(arhosek_tristim_skymodel_radiance(skyStates[i], theta, gamma, i));
                }

                // TODO: a better way to calculate luminance
                // approx luminous efficiency as 0.29 
                double luminousEfficiency = 0.29 * 683.0;

                RGBSpectrum luminance = radiance * float(luminousEfficiency);
                storePixel(p, t, luminance[0], luminance[1], luminance[2]);
            }
        }

        mOldEnvMap = mEnvMap;
        mEnvMap = Texture::create2D(nPhi, nTheta, ResourceFormat::RGBA32Float, 1, 1, pixelsData);

        for (int i = 0; i < _countof(skyStates); ++i)
        {
            arhosekskymodelstate_free(skyStates[i]);
        }

        delete[] pixelsData;
    }

    void SunLight::updateAsnyc(const InternalData data)
    {
        auto& func = [=]
        {
            updateLightInfo(data);
            updateEnvironmentMap(data);
        };

        static ThreadPool<16> sThreadPool;
        sThreadPool.getAvailable() = std::thread(func);
    }
}
