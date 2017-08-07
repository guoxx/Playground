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
#include "ToneMapping.h"
#include "API/RenderContext.h"
#include "Graphics/FboHelper.h"

namespace Falcor
{
    static const char* kShaderFilename = "Effects\\ToneMapping.ps.slang";
    const Gui::DropdownList kOperatorList = { 
    { (uint32_t)ToneMapping::Operator::Clamp, "Clamp to LDR" },
    { (uint32_t)ToneMapping::Operator::Linear, "Linear" }, 
    { (uint32_t)ToneMapping::Operator::Reinhard, "Reinhard" },
    { (uint32_t)ToneMapping::Operator::ReinhardModified, "Modified Reinhard" }, 
    { (uint32_t)ToneMapping::Operator::HejiHableAlu, "Heji's approximation" },
    { (uint32_t)ToneMapping::Operator::HableUc2, "Uncharted 2" },
    { (uint32_t)ToneMapping::Operator::Aces, "ACES" }
    };

    // @@guoxx
    const Gui::DropdownList kExposureModeList = {
        { (uint32_t)ExposureMode_Manual_SBS, "Saturation-based Speed" },
        { (uint32_t)ExposureMode_Manual_SOS, "Standard output sensitivity" },
        { (uint32_t)ExposureMode_Automatic, "Automatic" },
        { (uint32_t)ExposureMode_Falcor, "Falcor default" },
    };

    const Gui::DropdownList kShutterSpeedList = {
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over1, "1s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over2, "1/2s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over4, "1/4s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over8, "1/8s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over15, "1/15s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over30, "1/30s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over60, "1/60s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over125, "1/125s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over250, "1/250s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over500, "1/500s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over1000, "1/1000s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over2000, "1/2000s" },
        { (uint32_t)ToneMapping::ShutterSpeed::ShutterSpeed1Over4000, "1/4000s" },
    };

    const Gui::DropdownList kFStopList = {
        { (uint32_t)ToneMapping::FStop::FStop1Point8, "f/1.8" },
        { (uint32_t)ToneMapping::FStop::FStop2Point0, "f/2.0" },
        { (uint32_t)ToneMapping::FStop::FStop2Point2, "f/2.2" },
        { (uint32_t)ToneMapping::FStop::FStop2Point5, "f/2.5" },
        { (uint32_t)ToneMapping::FStop::FStop2Point8, "f/2.8" },
        { (uint32_t)ToneMapping::FStop::FStop3Point2, "f/3.2" },
        { (uint32_t)ToneMapping::FStop::FStop3Point5, "f/3.5" },
        { (uint32_t)ToneMapping::FStop::FStop4Point0, "f/4.0" },
        { (uint32_t)ToneMapping::FStop::FStop4Point5, "f/4.5" },
        { (uint32_t)ToneMapping::FStop::FStop5Point0, "f/5.0" },
        { (uint32_t)ToneMapping::FStop::FStop5Point6, "f/5.6" },
        { (uint32_t)ToneMapping::FStop::FStop6Point3, "f/6.3" },
        { (uint32_t)ToneMapping::FStop::FStop7Point1, "f/7.1" },
        { (uint32_t)ToneMapping::FStop::FStop8Point0, "f/8.0" },
        { (uint32_t)ToneMapping::FStop::FStop9Point0, "f/9.0" },
        { (uint32_t)ToneMapping::FStop::FStop10Point0, "f/10.0" },
        { (uint32_t)ToneMapping::FStop::FStop11Point0, "f/11.0" },
        { (uint32_t)ToneMapping::FStop::FStop13Point0, "f/13.0" },
        { (uint32_t)ToneMapping::FStop::FStop14Point0, "f/14.0" },
        { (uint32_t)ToneMapping::FStop::FStop16Point0, "f/16.0" },
        { (uint32_t)ToneMapping::FStop::FStop18Point0, "f/18.0" },
        { (uint32_t)ToneMapping::FStop::FStop20Point0, "f/20.0" },
        { (uint32_t)ToneMapping::FStop::FStop22Point0, "f/22.0" },
    };

    const Gui::DropdownList kISORatingList = {
        { (uint32_t)ToneMapping::ISORating::ISO100, "ISO100" },
        { (uint32_t)ToneMapping::ISORating::ISO200, "ISO200" },
        { (uint32_t)ToneMapping::ISORating::ISO400, "ISO400" },
        { (uint32_t)ToneMapping::ISORating::ISO800, "ISO800" },
    };

    inline float getApertureFNumber(int32_t v)
    {
        static const float FNumbers[] =
        {
            1.8f, 2.0f, 2.2f, 2.5f, 2.8f, 3.2f, 3.5f, 4.0f, 4.5f, 5.0f, 5.6f, 6.3f, 7.1f, 8.0f,
            9.0f, 10.0f, 11.0f, 13.0f, 14.0f, 16.0f, 18.0f, 20.0f, 22.0f
        };
        return FNumbers[int(v)];
    }

    inline float getShutterSpeedValue(int32_t v)
    {
        static const float ShutterSpeedValues[] =
        {
            1.0f / 1.0f, 1.0f / 2.0f, 1.0f / 4.0f, 1.0f / 8.0f, 1.0f / 15.0f, 1.0f / 30.0f,
            1.0f / 60.0f, 1.0f / 125.0f, 1.0f / 250.0f, 1.0f / 500.0f, 1.0f / 1000.0f, 1.0f / 2000.0f, 1.0f / 4000.0f,
        };
        return ShutterSpeedValues[int(v)];
    }

    inline float getISORatingValue(int32_t v)
    {
        static const float ISOValues[] =
        {
            100.0f, 200.0f, 400.0f, 800.0f
        };
        return ISOValues[int(v)];
    }

    ToneMapping::~ToneMapping() = default;

    ToneMapping::ToneMapping(ToneMapping::Operator op)
    {
        createLuminancePass();
        createToneMapPass(op);

        Sampler::Desc samplerDesc;
        samplerDesc.setFilterMode(Sampler::Filter::Point, Sampler::Filter::Point, Sampler::Filter::Point);
        mpPointSampler = Sampler::create(samplerDesc);
        samplerDesc.setFilterMode(Sampler::Filter::Linear, Sampler::Filter::Linear, Sampler::Filter::Point);
        mpLinearSampler = Sampler::create(samplerDesc);

        // @@guoxx
        // Sunny 16 rule : https://en.wikipedia.org/wiki/Sunny_16_rule
        mExposureMode = ExposureMode_Manual_SOS;
        mShutterSpeed = ShutterSpeed::ShutterSpeed1Over125;
        mAperture = FStop::FStop16Point0;
        mISO = ISORating::ISO100;
        mConstBufferData.camSettings.exposureMode = mExposureMode;
        mConstBufferData.camSettings.shutterSpeed = getShutterSpeedValue(int32_t(mShutterSpeed));
        mConstBufferData.camSettings.aperture = getApertureFNumber(int32_t(mAperture));
        mConstBufferData.camSettings.ISO = getISORatingValue(int32_t(mISO));
    }

    ToneMapping::UniquePtr ToneMapping::create(Operator op)
    {
        ToneMapping* pTM = new ToneMapping(op);
        return ToneMapping::UniquePtr(pTM);
    }

    void ToneMapping::createLuminanceFbo(Fbo::SharedPtr pSrcFbo)
    {
        bool createFbo = mpLuminanceFbo == nullptr;
        ResourceFormat srcFormat = pSrcFbo->getColorTexture(0)->getFormat();
        uint32_t bytesPerChannel = getFormatBytesPerBlock(srcFormat) / getFormatChannelCount(srcFormat);
        
        // Find the required texture size and format
        ResourceFormat luminanceFormat = (bytesPerChannel == 32) ? ResourceFormat::R32Float : ResourceFormat::R16Float;
        uint32_t requiredHeight = getLowerPowerOf2(pSrcFbo->getHeight());
        uint32_t requiredWidth = getLowerPowerOf2(pSrcFbo->getWidth());

        if(createFbo == false)
        {
            createFbo = (requiredWidth != mpLuminanceFbo->getWidth()) ||
                (requiredHeight != mpLuminanceFbo->getHeight()) ||
                (luminanceFormat != mpLuminanceFbo->getColorTexture(0)->getFormat());
        }

        if(createFbo)
        {
            Fbo::Desc desc;
            desc.setColorTarget(0, luminanceFormat);
            mpLuminanceFbo = FboHelper::create2D(requiredWidth, requiredHeight, desc, 1, Fbo::kAttachEntireMipLevel);
        }
    }

    void ToneMapping::execute(RenderContext* pRenderContext, Fbo::SharedPtr pSrc, Fbo::SharedPtr pDst)
    {
        GraphicsState::SharedPtr pState = pRenderContext->getGraphicsState();
        createLuminanceFbo(pSrc);

        //Set shared vars
        mpToneMapVars->setSrv(mBindLocations.colorTex.regSpace, mBindLocations.colorTex.baseRegIndex, 0, pSrc->getColorTexture(0)->getSRV());
        mpLuminanceVars->setSrv(mBindLocations.colorTex.regSpace, mBindLocations.colorTex.baseRegIndex, 0, pSrc->getColorTexture(0)->getSRV());
        mpToneMapVars->setSampler(mBindLocations.colorSampler.regSpace, mBindLocations.colorSampler.baseRegIndex, 0, mpPointSampler);
        mpLuminanceVars->setSampler(mBindLocations.colorSampler.regSpace, mBindLocations.colorSampler.baseRegIndex, 0, mpLinearSampler);

        //Calculate luminance
        pRenderContext->setGraphicsVars(mpLuminanceVars);
        pState->setFbo(mpLuminanceFbo);
        mpLuminancePass->execute(pRenderContext);
        mpLuminanceFbo->getColorTexture(0)->generateMips();

        //Set Tone map vars
        if (mOperator != Operator::Clamp)
        {
            mpToneMapCBuffer->setBlob(&mConstBufferData, 0u, sizeof(mConstBufferData));
            mpToneMapVars->setSampler(mBindLocations.luminanceSampler.regSpace, mBindLocations.luminanceSampler.baseRegIndex, 0, mpLinearSampler);
            mpToneMapVars->setSrv(mBindLocations.luminanceTex.regSpace, mBindLocations.luminanceTex.baseRegIndex, 0, mpLuminanceFbo->getColorTexture(0)->getSRV());
        }

        //Tone map
        pRenderContext->setGraphicsVars(mpToneMapVars);
        pState->setFbo(pDst);
        mpToneMapPass->execute(pRenderContext);
    }

    void ToneMapping::createToneMapPass(ToneMapping::Operator op)
    {
        mpToneMapPass = FullScreenPass::create(kShaderFilename);

        mOperator = op;
        switch(op)
        {
        case Operator::Clamp:
            mpToneMapPass->getProgram()->addDefine("_CLAMP");
            break;
        case Operator::Linear:
            mpToneMapPass->getProgram()->addDefine("_LINEAR");
            break;
        case Operator::Reinhard:
            mpToneMapPass->getProgram()->addDefine("_REINHARD");
            break;
        case Operator::ReinhardModified:
            mpToneMapPass->getProgram()->addDefine("_REINHARD_MOD");
            break;
        case Operator::HejiHableAlu:
            mpToneMapPass->getProgram()->addDefine("_HEJI_HABLE_ALU");
            break;
        case Operator::HableUc2:
            mpToneMapPass->getProgram()->addDefine("_HABLE_UC2");
            break;
        case Operator::Aces:
            mpToneMapPass->getProgram()->addDefine("_ACES");
            break;
        default:
            should_not_get_here();
        }

        const auto& pReflector = mpToneMapPass->getProgram()->getActiveVersion()->getReflector();
        mpToneMapVars = GraphicsVars::create(pReflector);
        mpToneMapCBuffer = mpToneMapVars["PerImageCB"];
        mBindLocations.luminanceSampler = getResourceBindLocation(pReflector.get(), "gLuminanceTexSampler");
        mBindLocations.colorSampler     = getResourceBindLocation(pReflector.get(), "gColorSampler");
        mBindLocations.colorTex         = getResourceBindLocation(pReflector.get(), "gColorTex");
        mBindLocations.luminanceTex     = getResourceBindLocation(pReflector.get(), "gLuminanceTex");

    }

    void ToneMapping::createLuminancePass()
    {
        mpLuminancePass = FullScreenPass::create(kShaderFilename);
        mpLuminancePass->getProgram()->addDefine("_LUMINANCE");
        const auto& pReflector = mpLuminancePass->getProgram()->getActiveVersion()->getReflector();
        mpLuminanceVars = GraphicsVars::create(pReflector);
    }

    void ToneMapping::renderUI(Gui* pGui, const char* uiGroup)
    {
        if((uiGroup == nullptr) || pGui->beginGroup(uiGroup))
        {
            uint32_t opIndex = static_cast<uint32_t>(mOperator);
            if (pGui->addDropdown("Operator", kOperatorList, opIndex))
            {
                mOperator = static_cast<Operator>(opIndex);
                createToneMapPass(mOperator);
            }

            pGui->addFloatVar("Exposure Key", mConstBufferData.exposureKey, 0.0001f, 200.0f);
            pGui->addFloatVar("Luminance LOD", mConstBufferData.luminanceLod, 0, 16, 0.025f);
            //Only give option to change these if the relevant operator is selected
            if (mOperator == Operator::ReinhardModified)
            {
                pGui->addFloatVar("White Luminance", mConstBufferData.whiteMaxLuminance, 0.1f, FLT_MAX, 0.2f);
            }
            else if (mOperator == Operator::HableUc2)
            {
                pGui->addFloatVar("Linear White", mConstBufferData.whiteScale, 0, 100, 0.01f);
            }

            // @@guoxx
            if (pGui->beginGroup("Camera Settings"))
            {
                if (pGui->addDropdown("Expsure Mode", kExposureModeList, mExposureMode))
                {
                    mConstBufferData.camSettings.exposureMode = mExposureMode;
                }

                uint32_t shutterSpeedIndex = static_cast<uint32_t>(mShutterSpeed);
                if (pGui->addDropdown("Shutter Speed", kShutterSpeedList, shutterSpeedIndex))
                {
                    mConstBufferData.camSettings.shutterSpeed = getShutterSpeedValue(shutterSpeedIndex);
                    mShutterSpeed = static_cast<ShutterSpeed>(shutterSpeedIndex);
                }

                uint32_t apertureIndex = static_cast<uint32_t>(mAperture);
                if (pGui->addDropdown("Aperture", kFStopList, apertureIndex))
                {
                    mConstBufferData.camSettings.aperture = getApertureFNumber(apertureIndex);
                    mAperture = static_cast<FStop>(apertureIndex);
                }

                uint32_t ISOIndex = static_cast<uint32_t>(mISO);
                if (pGui->addDropdown("ISO", kISORatingList, ISOIndex))
                {
                    mConstBufferData.camSettings.ISO = getISORatingValue(ISOIndex);
                    mISO = static_cast<ISORating>(ISOIndex);
                }

                pGui->endGroup();
            }

            if (uiGroup) pGui->endGroup();
        }
    }

    void ToneMapping::setOperator(Operator op)
    {
        if(op != mOperator)
        {
            createToneMapPass(op);
        }
    }

    void ToneMapping::setExposureKey(float exposureKey)
    {
        mConstBufferData.exposureKey = max(0.001f, exposureKey);
    }

    void ToneMapping::setWhiteMaxLuminance(float maxLuminance)
    {
        mConstBufferData.whiteMaxLuminance = maxLuminance;
    }

    void ToneMapping::setLuminanceLod(float lod)
    {
        mConstBufferData.luminanceLod = clamp(lod, 0.0f, 16.0f);
    }

    void ToneMapping::setWhiteScale(float whiteScale)
    {
        mConstBufferData.whiteScale = max(0.001f, whiteScale);
    }
}