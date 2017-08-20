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
__import Helpers;
__import ShaderCommon;

#include "CloudsCommon.hlsli"


Texture3D baseShapeLookup : register(t0);
Texture3D erosionLookup : register(t1);
Texture2D curlNoiseLookup : register(t2);
Texture2D weatherLookup : register(t3);

SamplerState texSampler : register(s0);


// ----
static const float cb_cloudSpeed = 10.0;
static const float cb_cloudTopOffset = 500.0;

static const float3 cb_windDirection = float3(1.0, 0.0, 0.0);

static const float3 cb_weatherTexMod = float3(16.0, 0, 0); // scale(x), offset(y, z)
static const float cb_windStrength = 0;
// ----

static float cloudBottomFade = 0.3;
static float highFreqNoiseScale = 8;
static float cloudDistortion = 0.45;
static float cloudDistortionScale = 0.5;
static const float cb_groundRadius = 6360000;
static const float cb_cloudVolumeStartHeight = 1500;
static const float cb_cloudVolumeHeight = 2500;

static const float VOLUME_END_HEIGHT = cb_cloudVolumeStartHeight + cb_cloudVolumeHeight;
// planet center (world space)
static const float3 PLANET_CENTER = float3(0.0f, -cb_groundRadius, 0.0f);
// radius from the planet center to the bottom of the cloud volume
static const float PLANET_CENTER_TO_LOWER_CLOUD_RADIUS = cb_groundRadius + cb_cloudVolumeStartHeight;
// radius from the planet center to the top of the cloud volume
static const float PLANET_CENTER_TO_UPPER_CLOUD_RADIUS = cb_groundRadius + VOLUME_END_HEIGHT;

static const float CLOUD_SCALE = 1.0f / VOLUME_END_HEIGHT;
static const float3 WEATHER_TEX_MOD = float3(1.0f / (VOLUME_END_HEIGHT * cb_weatherTexMod.x), cb_weatherTexMod.y, cb_weatherTexMod.z); 
static const float2 WEATHER_TEX_MOVE_SPEED = float2(cb_windStrength * cb_windDirection.x, cb_windStrength * cb_windDirection.z); // this is modded by app run time

// samples based on shell thickness between inner and outer volume
static const uint2 SAMPLE_RANGE = uint2(64u, 128u);

static const float4 STRATUS_GRADIENT = float4(0.011f, 0.098f, 0.126f, 0.225f);
static const float4 STRATOCUMULUS_GRADIENT = float4(0.0f, 0.096f, 0.311f, 0.506f);
static const float4 CUMULUS_GRADIENT = float4(0.0f, 0.087f, 0.749f, 1.0f);

// random vectors on the unit sphere
static const float3 RANDOM_VECTORS[] =
{
	float3( 0.38051305f,  0.92453449f, -0.02111345f),
	float3(-0.50625799f, -0.03590792f, -0.86163418f),
	float3(-0.32509218f, -0.94557439f,  0.01428793f),
	float3( 0.09026238f, -0.27376545f,  0.95755165f),
	float3( 0.28128598f,  0.42443639f, -0.86065785f),
	float3(-0.16852403f,  0.14748697f,  0.97460106f)
};

static const uint LIGHT_RAY_ITERATIONS = 6u;
static const float RCP_LIGHT_RAY_ITERATIONS = 1.0f / float(LIGHT_RAY_ITERATIONS);


/**
 * Perform a ray-sphere intersection test.
 * Returns the number of intersections in the direction of the ray (excludes intersections behind the ray origin), between 0 and 2.
 * In the case of more than one intersection, the nearest point will be returned in t1.
 * 
 * http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
 */
uint intersectRaySphere(
	float3 rayOrigin,
	float3 rayDir, // must be normalized
	float3 sphereCenter,
	float  sphereRadius,
	out float2 t)
{
	float3 l = rayOrigin - sphereCenter;
	float a = 1.0f; // dot(rayDir, rayDir) where rayDir is normalized
	float b = 2.0f * dot(rayDir, l);
	float c = dot(l, l) - sphereRadius * sphereRadius;
	float discriminate = b * b - 4.0f * a * c;
	if(discriminate < 0.0f)
	{
		t.x = t.y = 0.0f;
		return 0u;
	}
	else if(abs(discriminate) - 0.00005f <= 0.0f)
	{
		t.x = t.y = -0.5f * b / a;
		return 1u;
	}
	else
	{
		float q = b > 0.0f ?
			-0.5f * (b + sqrt(discriminate)) :
			-0.5f * (b - sqrt(discriminate));
		float h1 = q / a;
		float h2 = c / q;
		t.x = min(h1, h2);
		t.y = max(h1, h2);
		if(t.x < 0.0f)
		{
			t.x = t.y;
			if(t.x < 0.0f)
			{
				return 0u;
			}
			return 1u;
		}
		return 2u;
	}
}

float remap(float value, float oldMin, float oldMax, float newMin, float newMax)
{
	return newMin + (value - oldMin) / (oldMax - oldMin) * (newMax - newMin);
}

float3 sampleWeather(float3 pos)
{
#if 1
	float2 uv = pos.xz * WEATHER_TEX_MOD.x;// + WEATHER_TEX_MOD.yz + (WEATHER_TEX_MOVE_SPEED * gClouds.mGlobalTime);
	uv = uv * 0.5 + 0.5;
	float3 data = weatherLookup.SampleLevel(texSampler, uv, 0.0f).rgb;
	return float3(data.r, gClouds.mWeatherData.g, 1.0f - data.b);
#else
	float2 uv = pos.xz * WEATHER_TEX_MOD.x;// + WEATHER_TEX_MOD.yz + (WEATHER_TEX_MOVE_SPEED * gClouds.mGlobalTime);
	uv = uv * 0.5 + 0.5;
	return weatherLookup.SampleLevel(texSampler, uv, 0.0f).rgb;
#endif
}

float getCoverage(float3 weatherData)
{
	return weatherData.r;
}

float getPrecipitation(float3 weatherData)
{
	return weatherData.g;
}

float getCloudType(float3 weatherData)
{
	// weather b channel tells the cloud type 0.0 = stratus, 0.5 = stratocumulus, 1.0 = cumulus
	return weatherData.b;
}

float heightFraction(float3 pos)
{
	return saturate((distance(pos, PLANET_CENTER) - PLANET_CENTER_TO_LOWER_CLOUD_RADIUS) / cb_cloudVolumeHeight);
}

float gradientStep(float a, float4 gradient)
{
	return smoothstep(gradient.x, gradient.y, a) - smoothstep(gradient.z, gradient.w, a);
}

float4 mixGradients(float cloudType)
{
#if 1
	if (cloudType < 0.5)
	{
		return lerp(STRATUS_GRADIENT, STRATOCUMULUS_GRADIENT, cloudType * 2.0);
	}
	else
	{
		return lerp(STRATOCUMULUS_GRADIENT, CUMULUS_GRADIENT, (cloudType - 0.5) * 2.0);
	}
#else
	float stratus = 1.0f - saturate(cloudType * 2.0f);
	float stratocumulus = 1.0f - abs(cloudType - 0.5f) * 2.0f;
	float cumulus = saturate(cloudType - 0.5f) * 2.0f;
	return STRATUS_GRADIENT * stratus + STRATOCUMULUS_GRADIENT * stratocumulus + CUMULUS_GRADIENT * cumulus;
#endif
}

float densityHeightGradient(float heightFrac, float cloudType)
{
	float4 cloudGradient = mixGradients(cloudType);
#if 0
	return remap(heightFrac, cloudGradient.x, cloudGradient.y, 0, 1) * remap(heightFrac, cloudGradient.z, cloudGradient.w, 1, 0);
#else
	return smoothstep(cloudGradient.x, cloudGradient.y, heightFrac) - smoothstep(cloudGradient.z, cloudGradient.w, heightFrac);
#endif
}

float smoothThreshold(float value, float threshold, float edgeSize)
{
	return smoothstep(threshold, threshold + edgeSize, value);
}

float sampleCloudDensity(
	float3 pos,
	float3 weatherData,
	float heightFrac,
	float lod)
{
	// Skew in wind direction.
	pos += heightFrac * cb_windDirection * cb_cloudTopOffset;

	// Animate clouds in wind direction and add a small upward
	// bias to the wind direction
	pos += (cb_windDirection + float3(0.0f, 0.1f, 0.0f)) * cb_cloudSpeed * gClouds.mGlobalTime;

	// TODO: why?
	pos *= CLOUD_SCALE;

	// Read the low-frequency Perlin-Worley and Worley noises.	
	float4 lowFreqNoise = baseShapeLookup.SampleLevel(texSampler, pos, lerp(0.0f, gClouds.mBaseShapeTextureBottomMipLevel, lod));

    // Build an FBM out of the low frequency Worley noises
    // that can be used to add detail to the low-frequency
    // Perlin-Worley noise.
	float lowFreqFBM =
		(lowFreqNoise.g * 0.625f) +
		(lowFreqNoise.b * 0.25f) +
		(lowFreqNoise.a * 0.125f);

    // define the base cloud shape by dilating it with the
    // low-frequency FBM made of Worley noise.
	float baseCloud = remap(
		lowFreqNoise.r,
		-(1.0f - lowFreqFBM), 1.0f, // gets about the same results just using -lowFreqFBM
		0.0f, 1.0f);

    // Get the density-height gradient using the density height
    // function explained in Section 4.3.2.
	float densityGradient = densityHeightGradient(heightFrac, getCloudType(weatherData));

	// Apply the height function to the base cloud shape.
	baseCloud *= densityGradient;

	//baseCloud *= 1 - pow(heightFrac, 22);

	float cloudCoverage = getCoverage(weatherData);

	float anvil_bias = 0;
	cloudCoverage = pow(cloudCoverage, remap(1-heightFrac, 0.7, 0.8, 1.0, lerp(1.0, 0.5, anvil_bias)));

#if 0
	float baseCloudWithCoverage = saturate(remap(baseCloud, saturate(heightFrac / cloudCoverage), 1.0, 0.0, 1.0)) * cloudCoverage;
#else
	// Use remap to apply the cloud coverage attribute.
	float baseCloudWithCoverage = saturate(remap(
		baseCloud,
		1.0f - cloudCoverage, 1.0f,
		0.0f, 1.0f));

    // Multiply the result by the cloud coverage attribute so
    // that smaller clouds are lighter and more aesthetically
    // pleasing.
	baseCloudWithCoverage *= cloudCoverage;
#endif

	float finalCloud = baseCloudWithCoverage;

	if (gClouds.mEnableHighFreqNoise)
	{
		float3 coord = pos * highFreqNoiseScale;

		if (gClouds.mEnableCurlNoise)
		{
			// Add some turbulence to bottoms of clouds.
			float3 curlNoise = curlNoiseLookup.Sample(texSampler, pos.xy * cloudDistortionScale).rgb;
			curlNoise *= cloudDistortion * (1.0f - heightFrac);
			coord += curlNoise;
		}

		// Sample high-frequency noises.
		float3 highFreqNoise = erosionLookup.SampleLevel(texSampler, coord, lerp(0.0f, gClouds.mErosionTextureBottomMipLevel, lod)).gba;

		// Build high-frequency Worley noise FBM.
		float highFreqFBM =
			(highFreqNoise.r * 0.625f) +
			(highFreqNoise.g * 0.25f) +
			(highFreqNoise.b * 0.125f);

	    // Transition from wispy shapes to billowy shapes over height.
		float highFreqNoiseModifier = lerp(highFreqFBM, 1.0f - highFreqFBM, heightFrac);

	    // Erode the base cloud shape with the distorted
	    // high-frequency Worley noises.
		finalCloud = remap(
			baseCloudWithCoverage,
			highFreqNoiseModifier * 0.2f, 1.0f,
			0.0f, 1.0f);
	}

	if (gClouds.mEnableHeightFade)
	{
		finalCloud *= smoothstep(0.0, cloudBottomFade, heightFrac);
	}

	return finalCloud;
}

float beerLambert(float sampleDensity, float precipitation)
{
	return exp(-sampleDensity * precipitation);
}

float powder(float sampleDensity, float lightDotEye)
{
	float powd = 1.0f - exp(-sampleDensity * 2.0f);
	return lerp(
		1.0f,
		powd,
		saturate((-lightDotEye * 0.5f) + 0.5f) // [-1,1]->[0,1]
	);
}

float henyeyGreenstein(
	float lightDotEye,
	float g)
{
	float g2 = g * g;
	return ((1.0f - g2) / pow((1.0f + g2 - 2.0f * g * lightDotEye), 1.5f)) * 0.25f;
}

float lightEnergy(
	float lightDotEye,
	float densitySample,
	float originalDensity,
	float precipitation)
{
	return 2.0f *
		beerLambert(densitySample, precipitation) *
		powder(originalDensity, lightDotEye) * 
		lerp(henyeyGreenstein(lightDotEye, 0.8f), henyeyGreenstein(lightDotEye, -0.5f), 0.5f);
}

// TODO get from cb values - has to change as time of day changes
float3 ambientLight(float heightFrac)
{
	return lerp(
		float3(0.5f, 0.67f, 0.82f),
		float3(1.0f, 1.0f, 1.0f),
		heightFrac);
}

float sampleCloudDensityAlongCone(
	float3 startPos,
	float  stepSize,
	float  lightDotEye,
	float  originalDensity)
{
	float3 lightStep = stepSize * -gClouds.mSunLightDirection;
	float3 pos = startPos;
	float coneRadius = 1.0f;
	float coneStep = RCP_LIGHT_RAY_ITERATIONS;
	float densityAlongCone = 0.0f;
	float lod = 0.0f;
	float lodStride = RCP_LIGHT_RAY_ITERATIONS;
	float3 weatherData = 0.0f;
	float rcpThickness = 1.0f / (stepSize * LIGHT_RAY_ITERATIONS);
	float density = 0.0f;

	for(uint i = 0u; i < LIGHT_RAY_ITERATIONS; ++i)
	{
		float3 conePos = pos + coneRadius * RANDOM_VECTORS[i] * float(i + 1u);
		float heightFrac = heightFraction(conePos);
		if(heightFrac <= 1.0f)
		{
			weatherData = sampleWeather(conePos);
			float cloudDensity = sampleCloudDensity(
				conePos,
				weatherData,
				heightFrac,
				lod);
			if(cloudDensity > 0.0f)
			{
				density += cloudDensity;
				float transmittance = 1.0f - (density * rcpThickness);
				densityAlongCone += (cloudDensity * transmittance);
			}
		}
		pos += lightStep;
		coneRadius += coneStep;
		lod += lodStride;
	}
	// take additional step at large distance away for shadowing from other clouds
	pos = pos + (lightStep * 8.0f);
	weatherData = sampleWeather(pos);
	float heightFrac = heightFraction(pos);
	if(heightFrac <= 1.0f)
	{
		float cloudDensity = sampleCloudDensity(
			pos,
			weatherData,
			heightFrac,
			0.8f);
		// no need to branch here since density variable is no longer used after this
		density += cloudDensity;
		float transmittance = 1.0f - saturate(density * rcpThickness);
		densityAlongCone += (cloudDensity * transmittance);
	}

	return saturate(lightEnergy(
		lightDotEye,
		densityAlongCone,
		originalDensity,
		lerp(1.0f, 2.0f, getPrecipitation(weatherData))));
}

float4 traceClouds(
	float3 eyePosW,			// world space view position
	float3 viewDirW,		// world space view direction
	float3 startPos,		// world space start position
	float3 endPos)			// world space end position
{
	float3 dir = endPos - startPos;
	float thickness = length(dir);
	float rcpThickness = 1.0f / thickness;
	uint sampleCount = lerp(SAMPLE_RANGE.x, SAMPLE_RANGE.y, saturate((thickness - cb_cloudVolumeHeight) / cb_cloudVolumeHeight));
	float stepSize = thickness / float(sampleCount);
	dir /= thickness;
	float3 posStep = stepSize * dir;

	float lightDotEye = -dot(gClouds.mSunLightDirection, viewDirW);

	float3 pos = startPos;
	float3 weatherData = 0.0f;
	float4 result = 0.0f;
	float density = 0.0f;
	
	[loop]
	for(uint i = 0u; i < sampleCount; ++i)
	{
		float heightFrac = heightFraction(pos);
		weatherData = sampleWeather(pos);
		float cloudDensity = sampleCloudDensity(
			pos,
			weatherData,
			heightFrac,
			0.0f);

		if(cloudDensity > 0.0f)
		{
			density += cloudDensity;
			float transmittance = 1.0f - (density * rcpThickness);
			float lightDensity = sampleCloudDensityAlongCone(
				pos,
				stepSize,
				lightDotEye,
				cloudDensity);

			float3 ambientBadApprox = ambientLight(heightFrac) * min(1.0f, length(gClouds.mSunIrradiance.rgb * 0.0125f)) * transmittance;
			float4 source = float4((gClouds.mSunIrradiance.rgb * lightDensity) + ambientBadApprox, cloudDensity * transmittance); // TODO enable ambient when added to constant buffer
			source.rgb *= source.a;
			result = (1.0f - result.a) * source + result;
			//if(result.a >= 1.0f) break;
		}

		pos += posStep;
	}

#if 1
	return float4(result.xyz, 1 - saturate(result.a));
#else
	// experimental fog - may not be needed if clouds are drawn before atmosphere - would have to draw sun by itself, then clouds, then atmosphere
	// fogAmt = 0 to disable
	float fogAmt = 1.0f - exp(-distance(startPos, eyePosW) * 0.00001f);
	float3 fogColor = float3(0.3f, 0.4f, 0.45f) * length(gClouds.mSunIrradiance.rgb * 0.125f) * 0.8f;
	float3 sunColor = normalize(gClouds.mSunIrradiance.rgb) * 4.0f * length(gClouds.mSunIrradiance.rgb * 0.125f);
	fogColor = lerp(fogColor, sunColor, pow(saturate(lightDotEye), 8.0f));
	return float4(clamp(lerp(result.rgb, fogColor, fogAmt), 0.0f, 1000.0f), saturate(result.a));
#endif
}

float4 main(VsOut pIn) : SV_TARGET
{
	float4 posW = mul(float4(pIn.dir, 1), gClouds.mCamera.invViewProj);
	posW /= posW.w;

	float3 eyePosW = gClouds.mCamera.position;
	float3 viewDirW = normalize(posW.xyz - eyePosW);

	// @@guoxx: we are limited to ground view, so depth test is good enough for us
	bool depthPresent = false;

	// find nearest planet surface point
	float2 ph = 0.0f;
	uint planetHits = intersectRaySphere(
		eyePosW,
		viewDirW,
		PLANET_CENTER,
		cb_groundRadius,
		ph);

	// find nearest inner shell point
	float2 ih = 0.0f;
	uint innerShellHits = intersectRaySphere(
		eyePosW,
		viewDirW,
		PLANET_CENTER,
		PLANET_CENTER_TO_LOWER_CLOUD_RADIUS,
		ih);

	// find nearest outer shell point
	float2 oh = 0.0f;
	uint outerShellHits = intersectRaySphere(
		eyePosW,
		viewDirW,
		PLANET_CENTER,
		PLANET_CENTER_TO_UPPER_CLOUD_RADIUS,
		oh);

	// world space ray intersections
	float3 planetHitSpot = eyePosW + (viewDirW * ph.x);
	float3 innerShellHit = eyePosW + (viewDirW * ih.x);
	float3 outerShellHit = eyePosW + (viewDirW * oh.x);

	// eye radius from planet center
	float eyeRadius = distance(eyePosW, PLANET_CENTER);
	if(eyeRadius < PLANET_CENTER_TO_LOWER_CLOUD_RADIUS) // under inner shell
	{
		// exit if there's something in front of the start of the cloud volume
		if((depthPresent && (distance(posW.xyz, eyePosW) < distance(innerShellHit, eyePosW))) ||
			planetHits > 0u) // shell hits are guaranteed, but the ground may be occluding cloud layer
		{
			return float4(0.0f, 0.0f, 0.0f, 0.0f);
		}
		return traceClouds(
			eyePosW,
			viewDirW,
			innerShellHit,
			outerShellHit);
	}
	else if(eyeRadius > PLANET_CENTER_TO_UPPER_CLOUD_RADIUS) // over outer shell
	{
		// possibilities are
		// 1) enter outer shell, leave inner shell
		// 2) enter outer shell, leave outer shell
		float3 firstShellHit = outerShellHit;
		if(outerShellHits == 0u ||
			(depthPresent && (distance(posW.xyz, eyePosW) < distance(firstShellHit, eyePosW))))
		{
			return float4(0.0f, 0.0f, 0.0f, 0.0f);
		}
		float3 secondShellHit = outerShellHits == 2u && innerShellHits == 0u ? eyePosW + (viewDirW * oh.y) : innerShellHit;
		return traceClouds(
			eyePosW,
			viewDirW,
			firstShellHit,
			secondShellHit);
	}
	else // between shells
	{
		/*
		 * From a practical viewpoint (properly scaled planet, atmosphere, etc.)
		 * only one shell will be hit.
		 * Start position is always eye position.
		 */
		float3 shellHit = innerShellHits > 0u ? innerShellHit : outerShellHit;
		// if there's something in the depth buffer that's closer, that's the end point
		if(depthPresent && (distance(posW.xyz, eyePosW) < distance(shellHit, eyePosW)))
		{
			shellHit = posW.xyz;
		}
		return traceClouds(
			eyePosW,
			viewDirW,
			eyePosW,
			shellHit);
	}
}