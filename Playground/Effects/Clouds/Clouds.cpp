/***************************************************************************
# Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
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
#include "Clouds.h"
#include "glm/gtx/transform.hpp"
#include "Graphics/TextureHelper.h"
#include "Graphics/Camera/Camera.h"
#include "Graphics/SunLight.h"
#include "Utils/Gui.h"

namespace Falcor
{

    Clouds::UniquePtr Clouds::create()
    {
        UniquePtr pClouds = UniquePtr(new Clouds());
        if(pClouds->createResources() == false)
        {
            return nullptr;
        }
        return pClouds;
    }

    bool Clouds::createResources()
    {
        mpEffect = FullScreenPass::create("Effects\\Clouds.vs.slang", "Effects\\Clouds.ps.slang");
        mpProgram = std::dynamic_pointer_cast<GraphicsProgram, Program>(mpEffect->getProgram());
        mpVars = GraphicsVars::create(mpProgram->getActiveVersion()->getReflector());

        mpLowFreqNoisesTex = createTextureFromFile("Textures\\CloudsLowFrequencyNoises.dds", false, false);
        mpVars->setTexture("baseShapeLookup", mpLowFreqNoisesTex);

        mpHighFreqNoisesTex = createTextureFromFile("Textures\\CloudsHighFrequencyNoises.dds", false, false);
        mpVars->setTexture("erosionLookup", mpHighFreqNoisesTex);

        mpCurlNoisesTex = createTextureFromFile("Textures\\CurlNoise.png", false, false);
        mpVars->setTexture("curlNoiseLookup", mpCurlNoisesTex);

        mpWeatherTex = createTextureFromFile("Textures\\weather_data.png", false, false);
        mpVars->setTexture("weatherLookup", mpWeatherTex);

        Sampler::Desc samplerDesc;
        samplerDesc.setFilterMode(Sampler::Filter::Linear, Sampler::Filter::Linear, Sampler::Filter::Linear);
        samplerDesc.setMaxAnisotropy(8);
        Sampler::SharedPtr pSampler = Sampler::create(samplerDesc);
        mpVars->setSampler("texSampler", pSampler);

        DepthStencilState::Desc dsDesc;
        dsDesc.setDepthWriteMask(false).setDepthFunc(DepthStencilState::Func::LessEqual).setDepthTest(true);
        mpDsState = DepthStencilState::create(dsDesc);

        BlendState::Desc bsDesc;
        bsDesc.setRtBlend(0, true);
        bsDesc.setRtParams(0, BlendState::BlendOp::Add, BlendState::BlendOp::Add,
                           BlendState::BlendFunc::One, BlendState::BlendFunc::SrcAlpha,
                           BlendState::BlendFunc::Zero, BlendState::BlendFunc::One);
        for (uint32_t i = 1; i < Fbo::getMaxColorTargetCount(); ++i)
        {
            bsDesc.setRenderTargetWriteMask(i, false, false, false, false);
        }
        mpBlendState = BlendState::create(bsDesc);

        return true;
    }

    void Clouds::renderUI(Gui* pGui, const char* group)
    {
        if (!group || pGui->beginGroup(group))
        {
            pGui->addFloatVar("Coverage", mControls.mWeatherData.r, 0, 1);
            pGui->addFloatVar("Precipitation", mControls.mWeatherData.g, 0, 1);
            pGui->addFloatVar("CloudType", mControls.mWeatherData.b, 0, 1);

            pGui->addCheckBox("Enable High Frequency Noise", mControls.mEnableHighFreqNoise);
            pGui->addCheckBox("Enable Curl Noise", mControls.mEnableCurlNoise);
            pGui->addCheckBox("Enable Height Fade", mControls.mEnableHeightFade);

            if (group)
            {
                pGui->endGroup();
            }
        }
    }

    void Clouds::render(RenderContext* pRenderCtx, Camera* pCamera, SunLight* pSunLight, float globalTime)
    {
        pCamera->setIntoConstantBuffer(mpVars["PerFrameCB"].get(), "gClouds.mCamera");

        mpVars["PerFrameCB"]["gClouds.mGlobalTime"] = globalTime;
        mpVars["PerFrameCB"]["gClouds.mSunLightDirection"] = pSunLight->getWorldDirection();
        mpVars["PerFrameCB"]["gClouds.mSunIrradiance"] = pSunLight->getIntensity();;
        mpVars["PerFrameCB"]["gClouds.mBaseShapeTextureBottomMipLevel"] = mControls.mBaseShapeTextureBottomMipLevel;
        mpVars["PerFrameCB"]["gClouds.mErosionTextureBottomMipLevel"] = mControls.mErosionTextureBottomMipLevel;
        mpVars["PerFrameCB"]["gClouds.mWeatherData"] = mControls.mWeatherData;
        mpVars["PerFrameCB"]["gClouds.mSkyAmbient"] = pSunLight->GetSkyAmbient();
        mpVars["PerFrameCB"]["gClouds.mEnableHighFreqNoise"] = mControls.mEnableHighFreqNoise;
        mpVars["PerFrameCB"]["gClouds.mEnableCurlNoise"] = mControls.mEnableCurlNoise;
        mpVars["PerFrameCB"]["gClouds.mEnableHeightFade"] = mControls.mEnableHeightFade;

        pRenderCtx->pushGraphicsVars(mpVars);
        mpEffect->execute(pRenderCtx, mpDsState, mpBlendState);
        pRenderCtx->popGraphicsVars();
    }
}