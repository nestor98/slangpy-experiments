// SPDX-License-Identifier: Apache-2.0

import "slangpy";

float add_floats(float a, float b)
{
    return a + b;
}

float3 add_vectors(float3 a, float3 b)
{
    return a + b;
}

float4 sample_texture_at_uv(float2 uv, SamplerState sampler, Texture2D<float4> texture)
{
    return texture.SampleLevel(sampler, uv, 0);
}

/***********/

typealias vec2 = float2;
typealias vec3 = float3;
typealias vec4 = float4;

// Matrix types in GLSL are a more subtle issue, because they have different
// semantics from their HLSL/Slang equivalents in a few key ways:
//
// * The infix `*` operator always performs component-wise multiplication
//   in HLSL/Slang, but in GLSL it sometimes performs linear-algebraic
//   products (whenever we have matrix*matrix, vector*matrix, or matrix*vector).
//   HLSL/Slang require a distinct `mul()` function for those cases.
//
// * Because of differences in terminology and conventions, a linear-algebraic
//   product like `M*v` in GLSL is equivalent to `mul(v, M)` in HLSL/Slang
//   (note the reversed order of operands).
//
// * Constructing a matrix or vector from a single scalar consistently
//   replicates that scalar across all components/elements in HLSL/Slang,
//   but in GLSL it instead produces a diagonal matrix.
//
// These differences are not something we can surmount by defining the
// GLSL matrix types as aliases of the standard Slang ones, so instead we
// must define the GLSL matrix types as wrappers around the Slang ones.
//
struct mat2
{
    float2x2 data;

    __init(float e00, float e01, float e10, float e11)
    {
        data = float2x2(e00, e01, e10, e11);
    }

    // TODO: We need to fill in the other intializers and members
    // available on matrices here.
};

void mulAssign(inout vec2 left, mat2 right)
{
    left = mul(right.data, left);
}

float fract(float value)
{
    return frac(value);
}

// a * b in GLSL!
float2 mul(float2 b, mat2 a) {
    return mul(b, a.data);
}

float2 operator*(mat2 a, float2 b)
{
    return mul(b, a); // have to be flipped!
}

struct mat3 {
    float3x3 data;

    __init(float e00, float e01, float e02,
           float e10, float e11, float e12,
           float e20, float e21, float e22)
    {
        data = float3x3(e00, e01, e02,
                        e10, e11, e12,
                        e20, e21, e22);
    }

    __init(vec3 a, vec3 b, vec3 c) {

        data = (float3x3(a, b, c));
    }
    // TODO: We need to fill in the other intializers and members
    // available on matrices here.
};

// a * b in GLSL!
float3 mul(float3 b, mat3 a) {
    return mul(b, a.data);
}

float3 operator *(mat3 a, float3 b)
{
    return mul(b, a); // have to be flipped!
}

