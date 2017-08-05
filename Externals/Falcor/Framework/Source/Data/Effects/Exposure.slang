#ifndef _FALCOR_EXPOSURE_H_
#define _FALCOR_EXPOSURE_H_

#include "HostDeviceData.h"

// The two functions below were based on code and explanations provided by Padraic Hennessy (@PadraicHennessy).
// See this for more info: https://placeholderart.wordpress.com/2014/11/21/implementing-a-physically-based-camera-manual-exposure/

/*
* Get an exposure using the Saturation-based Speed method.
*/
float getSaturationBasedExposure(float aperture,
                                 float shutterSpeed,
                                 float iso)
{
    float l_max = (7800.0f / 65.0f) * (aperture * aperture) / (iso * shutterSpeed);
    return 1.0f / l_max;
}
 
/*
* Get an exposure using the Standard Output Sensitivity method.
* Accepts an additional parameter of the target middle grey.
*/
float getStandardOutputBasedExposure(float aperture,
                                     float shutterSpeed,
                                     float iso,
                                     float middleGrey = 0.18f)
{
    float l_avg = (1000.0f / 65.0f) * (aperture * aperture) / (iso * shutterSpeed);
    return middleGrey / l_avg;
}

float computeEV100FromAvgLuminance(float avgLuminance)
{
    // We later use the middle gray at 12.7% in order to have
    // a middle gray at 18% with a sqrt (2) room for specular highlights
    // But here we deal with the spot meter measuring the middle gray
    // which is fixed at 12.5 for matching standard camera
    // constructor settings (i.e. calibration constant K = 12.5)
    // Reference : http://en.wikipedia.org/wiki/Film_speed
    return log2(avgLuminance * 100.0f / 12.5f);
}

float convertEV100ToExposure(float EV100)
{
    // Compute the maximum luminance possible with H_sbs sensitivity
    // maxLum = 78 / ( S * q ) * N^2 / t
    // = 78 / ( S * q ) * 2^ EV_100
    // = 78 / (100 * 0.65) * 2^ EV_100
    // = 1.2 * 2^ EV
    // Reference : http://en.wikipedia.org/wiki/Film_speed
    float maxLuminance = 1.2f * pow(2.0f, EV100);
    return 1.0f / maxLuminance;
}

float log2Exposure(CameraSettings cameraSettings, float avgLuminance)
{
    float exposure = 0.0f;

    if (cameraSettings.exposureMode == ExposureMode_Manual_SBS)
    {
        exposure = log2(getSaturationBasedExposure(cameraSettings.aperture, cameraSettings.shutterSpeed, cameraSettings.ISO));
    }
    else if (cameraSettings.exposureMode == ExposureMode_Manual_SOS)
    {
		float middleGrey = 0.18f;
        exposure = log2(getStandardOutputBasedExposure(cameraSettings.aperture, cameraSettings.shutterSpeed, cameraSettings.ISO, middleGrey));
    }
    else
    {
        float autoEV100 = computeEV100FromAvgLuminance(avgLuminance);
        exposure = log2(convertEV100ToExposure(autoEV100));
    }

    return exposure;
}

// Determines the color based on exposure settings
float3 getExposedColor(CameraSettings cameraSettings, float3 color, float avgLuminance, float offset, out float exposure)
{
    exposure = log2Exposure(cameraSettings, avgLuminance);
    exposure += offset;
    return exp2(exposure) * color;
}

#endif	// _FALCOR_EXPOSURE_H_