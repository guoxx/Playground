#pragma once
#include "glm/detail/type_vec3.hpp"

namespace Falcor
{
    class SphericalCoordinates
    {
    public:
        // Right hand coordinate
        // (1, 0, 0) -> North
        // (0, 1, 0) -> Zenith
        // (0, 0, 1) -> East
        static SphericalCoordinates FromSphere(glm::vec3 coord);
        static SphericalCoordinates FromThetaAndPhi(float theta, float phi);
        static glm::vec3 ToSphere(SphericalCoordinates coord);

        float GetZenith() const;

        float GetElevation() const;

        float GetAzimuth() const;

    private:
        float m_Elevation{ 0 };
        float m_Azimuth{ 0 };
    };
}
