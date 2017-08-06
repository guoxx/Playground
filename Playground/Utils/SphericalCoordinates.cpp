#include "SphericalCoordinates.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "glm/glm.hpp"

using namespace Falcor;

SphericalCoordinates SphericalCoordinates::FromSphere(glm::vec3 coord)
{
    glm::vec3 normalizedCoord = glm::normalize(coord);

    // [0, PI]
    float zenith = std::acos(normalizedCoord.y);

    // [-PI/2, PI/2]
    float elevation = float(M_PI_2 - zenith);

    // [0, 2*PI]
    float azimuth = std::atan2(normalizedCoord.z, normalizedCoord.x);
    if (azimuth < 0)
    {
        azimuth += float(M_PI * 2);
    }

    SphericalCoordinates sphericalCoord;
    sphericalCoord.m_Elevation = elevation;
    sphericalCoord.m_Azimuth= azimuth;
    return sphericalCoord;
}

SphericalCoordinates SphericalCoordinates::FromThetaAndPhi(float theta, float phi)
{
    assert(0 <= theta && theta <= M_PI);
    assert(0 <= phi && phi <= M_PI * 2);

    SphericalCoordinates sphericalCoord;
    sphericalCoord.m_Elevation = float(M_PI_2 - theta);
    sphericalCoord.m_Azimuth = phi;
    return sphericalCoord;
}

glm::vec3 SphericalCoordinates::ToSphere(SphericalCoordinates coord)
{
    float y = std::cos(coord.GetZenith());
    float x = std::sin(coord.GetZenith()) * std::cos(coord.GetAzimuth());
    float z = std::sin(coord.GetZenith()) * std::sin(coord.GetAzimuth());
    return glm::vec3(x, y, z);
}

float SphericalCoordinates::GetZenith() const
{
    return float(M_PI_2 - m_Elevation);
}

float SphericalCoordinates::GetElevation() const
{
    return m_Elevation;
}

float SphericalCoordinates::GetAzimuth() const
{
    return m_Azimuth;
}
