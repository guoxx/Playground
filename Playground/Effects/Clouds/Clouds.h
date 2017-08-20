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
#pragma once
#include <memory>
#include "API/Sampler.h"
#include "API/Texture.h"
#include "graphics/Model/Model.h"
#include "Graphics/Program.h"
#include "API/ConstantBuffer.h"
#include "API/DepthStencilState.h"
#include "API/RasterizerState.h"
#include "API/BlendState.h"

namespace Falcor
{
    class Gui;
    class SunLight;
    class RenderContext;

    class CloudsControl
    {
    public:
        uint32_t mBaseShapeTextureBottomMipLevel = 8u;
        uint32_t mErosionTextureBottomMipLevel = 8u;

        bool mEnableHighFreqNoise = true;
        bool mEnableCurlNoise = true;
        bool mEnableHeightFade = true;

        glm::vec3 mWeatherData = glm::zero<vec3>();
    };

    class Clouds
    {
    public:
        using UniquePtr = std::unique_ptr<Clouds>;

        static UniquePtr create();

        void Clouds::renderUI(Gui* pGui, const char* group);

        void render(RenderContext* pRenderCtx, Camera* pCamera, SunLight* pSunLight, float globalTime);

    private:
        Clouds() = default;
        bool createResources();

        Texture::SharedPtr mpLowFreqNoisesTex;
        Texture::SharedPtr mpHighFreqNoisesTex;
        Texture::SharedPtr mpCurlNoisesTex;
        Texture::SharedPtr mpWeatherTex;

        FullScreenPass::UniquePtr mpEffect;
        GraphicsProgram::SharedPtr mpProgram;
        GraphicsVars::SharedPtr mpVars;

        BlendState::SharedPtr mpBlendState;
        DepthStencilState::SharedPtr mpDsState;

        CloudsControl mControls;
    };
}