// from https://github.com/NVIDIAGameWorks/Falcor/blob/eb540f6748774680ce0039aaf3ac9279266ec521/Source/RenderPasses/SimplePostFX/SimplePostFX.cs.slang

/***************************************************************************
 # Copyright (c) 2015-23, NVIDIA CORPORATION. All rights reserved.
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
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
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
 **************************************************************************/

/**
 * Simple Post FX
 *
 * applies vignette, chromatic aberration, barrel distortion, brightness-dependant (de)saturation, offset/scale/power tint per pixel
 *
 * also applies bloom which is an energy preserving weighted sum of 'gaussians' plus an optional star shape.
 * the bloom is calculated by creating an image pyramid fine to coarse, using a 10x10 binomial kernel (implemented via 5x5 bilinear taps)
 * the downsampling kernel needs to be quite wide to avoid aliasing artefacts
 * after 8 downsampling steps (so that the final PSF is quite wide - at least 64*10 pixels diameter), we upsample coarse to fine,
 * at each level blending a small fraction (the 'bloom amount') of the coarser layer with the 'current' layer.
 * in this way, the final kernel is a blend of gaussians with doubling variance, plus a spike at the middle.
 * optionally, 6 star like lobes are added by brute force line sampling.
 */

#include "Utils/Color/ColorHelpers.slang"

// cbuffer PerFrameCB
struct PerFrameCB{
    uint2 gResolution;
    float2 gInvRes;
    float gVignetteAmount;
    float gChromaticAberrationAmount;
    float2 gBarrelDistort;
    float3 gSaturationCurve;
    float3 gColorOffset;
    float3 gColorScale;
    float3 gColorPower;
    float gWipe;
    bool gTonemap;
}

struct BloomCB {
    uint2 gResolution;
    float2 gInvRes;
    float gStar;
    float2 gStarDir1;
    float2 gStarDir2;
    float2 gStarDir3;
    float gBloomThreshold;
    float gBloomStrength;
    float gBloomClamp;
    float gBloomAmount;
    bool gInPlace;   
}

// https://github.com/shader-slang/slangpy-samples/blob/ed2bea680eb2388639ea2e0183e5c3f26d3b5162/experiments/mipmap/app.slang#L27
float3 tonemap_aces_film(float3 input)
{
    float3 x = input.xyz;
    float a = 2.51;
    float b = 0.03;
    float c = 2.43;
    float d = 0.59;
    float e = 0.14;
    float3 col = saturate((x * (a * x + b)) / (x * (c * x + d) + e));
    return col;
}

// Texture2D<float4> gSrc;
// Texture2D<float4> gBloomed;
// RWTexture2D<float4> gDst;
// SamplerState gSampler;

float4 blurFilter2x2(Texture2D<float4> src, SamplerState gSampler, const float2 uv, float uscale, float vscale)
{
    // 2x2 bilinear-tap kernel, effective footprint 4x4, with offsets chosen to sample from the binomial 1 3 3 1
    uscale *= 0.75f;
    vscale *= 0.75f;
    // clang-format off
    return 0.25f * (
        src.SampleLevel(gSampler, uv + float2(-uscale, -vscale), 0.f) +
        src.SampleLevel(gSampler, uv + float2(uscale, -vscale), 0.f) +
        src.SampleLevel(gSampler, uv + float2(-uscale, vscale), 0.f) +
        src.SampleLevel(gSampler, uv + float2(uscale, vscale), 0.f));
    // clang-format on
}

float4 blurFilter3x3(Texture2D<float4> src, SamplerState gSampler, const float2 uv, const float uscale, const float vscale)
{
    // 3x3 bilinear-tap kernel, effective footprint 6x6, with offsets chosen to sample from the binomial 1 5 10 10 5 1
    const float w0 = (10.f + 10.f) / 32.f;
    const float x1 = 1.5f + 1.f / (1.f + 5.f), w1 = (1.f + 5.f) / 32.f;
    const float x1u = x1 * uscale;
    const float x1v = x1 * vscale;
    // clang-format off
    return
        (src.SampleLevel(gSampler, uv + float2(-x1u, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(-x1u, x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, x1v), 0.f)) * (w1 * w1) +
        (src.SampleLevel(gSampler, uv + float2(0.f, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(0.f, x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(-x1u, 0.f), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, 0.f), 0.f)) * (w1 * w0) +
        (src.SampleLevel(gSampler, uv + float2(0.f, 0.f), 0.f)) * (w0 * w0);
    // clang-format on
}

float4 blurFilter5x5(Texture2D<float4> src, SamplerState gSampler, const float2 uv, const float uscale, const float vscale)
{
    // 5x5 bilinear-tap kernel, effective footprint 10x10, with offsets chosen to sample from the binomial 1 9 36 84 126 126 84 32 9 1
    const float w0 = (126.f + 126.f) / 512.f;
    const float x1 = 1.5f + 36.f / (36.f + 84.f), w1 = (36.f + 84.f) / 512.f;
    const float x2 = 3.5f + 1.f / (1.f + 9.f), w2 = (1.f + 10.f) / 512.f;
    const float x1u = x1 * uscale, x2u = x2 * uscale;
    const float x1v = x1 * vscale, x2v = x2 * vscale;
    // clang-format off
    float4 s0 =
        (src.SampleLevel(gSampler, uv + float2(-x2u, -x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x2u, -x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(-x2u, x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x2u, x2v), 0.f)) * (w2 * w2) +
        (src.SampleLevel(gSampler, uv + float2(-x1u, -x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, -x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(-x1u, x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, x2v), 0.f)) * (w1 * w2) +
        (src.SampleLevel(gSampler, uv + float2(0.f, -x2v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(0.f, x2v), 0.f)) * (w0 * w2);
    float4 s1 =
        (src.SampleLevel(gSampler, uv + float2(-x2u, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x2u, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(-x2u, x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x2u, x1v), 0.f)) * (w2 * w1) +
        (src.SampleLevel(gSampler, uv + float2(-x1u, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(-x1u, x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, x1v), 0.f)) * (w1 * w1) +
        (src.SampleLevel(gSampler, uv + float2(0.f, -x1v), 0.f) +
         src.SampleLevel(gSampler, uv + float2(0.f, x1v), 0.f)) * (w0 * w1);
    float4 s2 =
        (src.SampleLevel(gSampler, uv + float2(-x2u, 0.f), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x2u, 0.f), 0.f)) * (w2 * w0) +
        (src.SampleLevel(gSampler, uv + float2(-x1u, 0.f), 0.f) +
         src.SampleLevel(gSampler, uv + float2(x1u, 0.f), 0.f)) * (w1 * w0) +
        (src.SampleLevel(gSampler, uv + float2(0.f, 0.f), 0.f)) * (w0 * w0);
    // clang-format on
    return s0 + s1 + s2;
}

// apply a threshold to the color, returning black if below threshold, otherwise scales 
float4 thresholdBloom(float4 c, float threshold)
{
    float v = saturate((luminance(c.rgb) - threshold) / max(1e-4f, threshold));
    return float4(c.rgb * v, 1.);
}

void downsample(uint2 pixel, BloomCB cb, Texture2D<float4> gSrc, 
                RWTexture2D<float4> gDst, SamplerState gSampler)
{
    if (any(pixel.xy >= cb.gResolution))
        return;
    float2 uv = (pixel + 0.5f) * cb.gInvRes;
    float4 c = blurFilter5x5(gSrc, gSampler, uv, cb.gInvRes.x * 0.5f, cb.gInvRes.y * 0.5f);
    gDst[pixel] = thresholdBloom(c, cb.gBloomThreshold);
}

// [numthreads(16, 16, 1)]
void upsample(uint2 pixel, BloomCB cb,
              Texture2D<float4> gSrc,
              Texture2D<float4> gBloomed,
              RWTexture2D<float4> gDst,
              SamplerState gSampler)
{
    if (any(pixel.xy >= cb.gResolution))
        return;
    float2 uv = (pixel + 0.5f) * cb.gInvRes;
    float4 centertap = blurFilter3x3(gBloomed, gSampler, uv, cb.gInvRes.x, cb.gInvRes.y);

    centertap.rgb = min(centertap.rgb * cb.gBloomStrength, cb.gBloomClamp); // scale up the bloom, and clamp to avoid excessive values
    // do a smooth clamping again instead:
    // centertap = thresholdBloom(centertap * cb.gBloomStrength, cb.gBloomClamp);
    if (cb.gStar > 0.f)
    {
        float4 star = 0.f;
        for (int i = -20; i <= 20; ++i)
        {
            float bri = 1.f / (10 + i * i);
            float4 c = gBloomed.SampleLevel(gSampler, uv + cb.gStarDir1 * i * 2, 0.f) +
                       gBloomed.SampleLevel(gSampler, uv + cb.gStarDir2 * i * 2, 0.f) +
                       gBloomed.SampleLevel(gSampler, uv + cb.gStarDir3 * i * 2, 0.f);
            star += c * bri;
        }
        star *= 1.f / 1.793358534742959f; // normalize by the total of the weights of star taps
        gDst[pixel] = lerp(lerp(gDst[pixel], centertap, cb.gBloomAmount), star, cb.gStar);
    }
    else
    {
        if (cb.gInPlace)
        {
            gDst[pixel] = lerp(gDst[pixel], centertap, cb.gBloomAmount);
        }
        else
        {
            gDst[pixel] = lerp(gSrc[pixel], centertap, cb.gBloomAmount);
        }
    }
}


void runPostFX(int2 pixel,
                 PerFrameCB cb,
                 Texture2D<float4> gSrc,
                 Texture2D<float4> gBloomed,
                 RWTexture2D<float4> gDst,
                 SamplerState gSampler)
{
    // float v = pixel.x / float(cb.gResolution.x);
    // gDst[pixel] = float4(v, .3, .6, 1.);
    if (any(pixel.xy >= cb.gResolution))
        return;
    const uint2 pixel = pixel.xy;
    float2 uv_circular = ((pixel + 0.5f) * 2.f - cb.gResolution.xy) * cb.gInvRes.y; // -1 to 1 vertically; more than -1 to 1 horizontally.
    float r_squared = dot(uv_circular, uv_circular);

    float barrel = (1.f + r_squared * cb.gBarrelDistort.y) * cb.gBarrelDistort.x;
    uv_circular *= barrel;

    float3 col;
    float2 ndc = float2(cb.gResolution.y * cb.gInvRes.x * 0.5f, 0.5f);
    float ca = 1.f + cb.gChromaticAberrationAmount;
    if (ca != 1.f)
    {
        col.x = gBloomed.SampleLevel(gSampler, uv_circular * ndc / ca + 0.5f, 0.f).x;
        col.y = gBloomed.SampleLevel(gSampler, uv_circular * ndc + 0.5f, 0.f).y;
        col.z = gBloomed.SampleLevel(gSampler, uv_circular * ndc * ca + 0.5f, 0.f).z;
    }
    else
    {
        col = gBloomed.SampleLevel(gSampler, uv_circular * ndc + 0.5f, 0.f).xyz;
    }
    // eyeballing the illuminance graphs with radius in https://grail.cs.washington.edu/projects/vignette/vign.iccv05.pdf,
    // it seems like there's quite a variety of falloffs with radius. let's go with a cheap circular falloff for now.
    float vignette = max(0.f, 1.f - r_squared * cb.gVignetteAmount);
    col *= vignette;

    float luma = luminance(col);
    // squashed_luma is used to partition the luminance range into a shadow-mid-bright, which we then fit a quadratic curve thru the desired
    // saturation
    float squashed_luma = saturate(luma / (luma + 0.5f));
    float saturation =
        max(0.f, cb.gSaturationCurve.z + squashed_luma * cb.gSaturationCurve.y + squashed_luma * squashed_luma * cb.gSaturationCurve.x);
    col = lerp(luma, col, saturation);

    // ASC-CDL style offset,power,slope (as implemented in blender, see
    // https://blender.stackexchange.com/questions/55231/what-is-the-the-asc-cdl-node)
    col = pow(max(0.f, col * cb.gColorScale + cb.gColorOffset), cb.gColorPower);
    
    float4 output = float4(col, 1.);
    if (pixel.x < cb.gWipe)
    {
        // output = gSrc.SampleLevel(gSampler, float2(pixel) * cb.gInvRes, 0.);
        output = gSrc[pixel];
    }

    if (cb.gTonemap) {
        output.xyz = tonemap_aces_film(output.xyz);
    }

    //output = float4(v, .3 * cb.gWipe * cb.gInvRes.x, .6, 1.);//gSrc[pixel];
    gDst[pixel] = output;
}