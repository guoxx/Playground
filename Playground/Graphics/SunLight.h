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
#pragma once
#include <string>
#include <glm/common.hpp>
#include "glm/geometric.hpp"
#include "API/Texture.h"
#include "glm/mat4x4.hpp"
#include "Data/HostDeviceData.h"
#include "Utils/Gui.h"
#include "Graphics/Model/Model.h"
#include "Graphics/Paths/MovableObject.h"
#include "Graphics/Light.h"
#include "Utils/SphericalCoordinates.h"
#include "Utils/Spectrum.h"

namespace Falcor
{
    class ConstantBuffer;
    class Gui;

    /** Directional light source.
    */
    class SunLight : public DirectionalLight, public std::enable_shared_from_this<SunLight>
    {
    public:
        using SharedPtr = std::shared_ptr<SunLight>;
        using SharedConstPtr = std::shared_ptr<const SunLight>;

        static SharedPtr create();

        SunLight();
        ~SunLight();

        /** Render UI elements for this light.
        \param[in] pGui The GUI to create the elements with
        \param[in] group Optional. If specified, creates a UI group to display elements within
        */
        void renderUI(Gui* pGui, const char* group = nullptr) override;

        /** Prepare GPU data
        */
        void prepareGPUData() override;

        /** Unload GPU data
        */
        void unloadGPUData() override;

        /** IMovableObject interface
        */
        void move(const glm::vec3& position, const glm::vec3& target, const glm::vec3& up) override;

        Texture::SharedPtr GetSkyEnvMap() const;

        glm::vec3 GetSkyAmbient() const;

    private:
        static SampledSpectrum computeSunRadiance(float sunTheta, float sunPhi, float turbidity, glm::vec3 groundAlbedo);

        struct InternalData
        {
            float mTurbidity;
            float mTheta;
            float mPhi;
            vec3 mGroundAlbedo; 
        } mData;

        void updateLightInfo(const InternalData& data);

        void updateEnvironmentMap(const InternalData& data);

        void updateAsnyc(const InternalData data);

        Texture::SharedPtr mEnvMap;
        glm::vec3 mSkyAmbient;
    };
}
