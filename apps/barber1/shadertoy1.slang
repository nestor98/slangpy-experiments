// https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang
// ^ Adapted this removing the interface, which was a lot of bloat for not a lot of functionality
// This file contains the pre definitions for glsl-slang compatibility
// the final slang will be a sandwich with shadertoy1.slang + user code (beginning of shadertoy2) + shadertoy2.slang

// some missing defines
static const int HW_PERFORMANCE = 1;

// ----- original header comments -----


// At this point we've described how our module will
// execute shader toy effects that implement the
// required interface, but we also need to set up
// the services that those effects are able to use.
//
// The shader toy "How To" file describes a large number of uniform
// shader parameters that are implicitly visible to every effect.
//
// If we were able to design an interface from scratch, we might
// prefer to make the `mainImage` function take some kind of
// explicit context parameter that provides access to these
// values, but because our goal is to be compatible with existing
// effects with their established `mainImage` signature, we will
// instead define these parameters using old-fashioned global-scope
// shader parameters.
//
cbuffer ShaderToyUniforms
{
    // Note: We do not currently define all of the parameters
    // exposed by Shader Toy, but rather just the most commonly
    // used ones.
    //
    // TODO: We can and should fill in the rest over time.
    //
    float4 iMouse, 
           iMouseRight, // right button mouse position
           iMouseShift, // some additional mouse positions of last time these keys were down
           iMouseCtrl,
           iMouseAlt;
    float2 iResolution;
    float iTime;

    float g_twist_mult;
    float g_freq_mult;
};

// In addition to the above parameters that use ordinary data types,
// shader toy also exposes the `iChannel*` parameters (`iChannel0`
// through `iChannel4`). These parameters represent sampled image
// inputs that can be bound to selected images as part of an effect.
//
// Traditional GLSL "sampler" types include both the texture image
// and sampler state, while Slang (like D3D, Vulkan, etc.) has
// distinct texture and sampler types. In order to define the
// channel variables in a way that is compatible with shader toy,
// we will define a `struct` type for a pair of a texture and
// a sampler:
//
struct TextureSamplerPair
{
    Texture2D t;
    SamplerState s;
};

// With our texture-sampler pair type defined, we can introduce
// the variables for the texture channels easily.
//
TextureSamplerPair iChannel0;
TextureSamplerPair iChannel1;
TextureSamplerPair iChannel2;
TextureSamplerPair iChannel3;

// TODO: Shader toy supports more than just 2D textures, so a good
// avenue for extension of the module would be to define an interface
// for the texture channels, and have implementations using various
// forms of textures.
//
// A really ambitious idea would be include one example of the
// texture-channel interface that uses an existing ShaderToy as a
// procedural texture.

// Shader toy effects access the contents of the `iChannel*` variables
// using the `texture()` function, so we need to provide a definition
// that is suitable:
//
float4 texture(TextureSamplerPair p, float3 uvw)
{
    // TODO: The right implementation to use here (at least
    // in the context of fragment shaders) is:
    //
    //      return p.t.Sample(p.s, uvw.xy);
    //
    // However, the current implementation of the main
    // application code doesn't include texture image
    // loading, so we will instead just fill in a
    // placeholder result for "texture" lookup:
    //
    return 0.5;
}

float4 texture(TextureSamplerPair p, float2 uv)
{
    return texture(p, float3(uv, 0));
}
// The last major issue we need to address in this module is the way that
// shader toy effects are authored in GLSL, which has several differences
// from Slang that could cause problems.
//
// Some of these differences can be surmounted relatively easily. For
// example, GLSL uses different names for its built-in vector types, but
// for the most part they are compatible with those defined by HLSL/Slang.
// We can paper over this difference by defining a few helpful type
// aliases.
//
typealias vec2 = float2;
typealias vec3 = float3;
typealias vec4 = float4;

typealias ivec2 = int2;
typealias ivec3 = int3;
typealias ivec4 = int4;
typealias bvec2 = bool2;
typealias bvec3 = bool3;

vec2 mod(vec2 x, float y)
{
    return fmod(x, y);
}

vec3 mod(vec3 x, float y)
{
    return fmod(x, y);
}

float mod(float x, float y)
{
    return fmod(x, y);
}

vec2 mix(vec2 x, vec2 y, float a)
{
    return lerp(x, y, a);
}
vec2 mix(vec2 x, vec2 y, vec2 a)
{
    return lerp(x, y, a);
}
vec3 mix(vec3 x, vec3 y, float a)
{
    return lerp(x, y, a);
}
vec4 mix(vec4 x, vec4 y, float a)
{
    return lerp(x, y, a);
}

float mix(float x, float y, float a)
{
    return lerp(x, y, a);
}

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

// TODO: fill in `mat3` and `mat4`.

// TODO: Ideally we would want to define overloaded operation functions
// to allow `*`, `*=`, etc. to apply to our user-space matrix types.
//
// Unfortunately, implementation bugs in the Slang compiler mean that
// user-defined operator overloads aren't working right now.
//
//      vec2 operator*(vec2 left, mat2 right)
//      {
//          return mul(right.data, left);
//      }
//
//      void operator*=(inout vec2 left, mat2 right)
//      {
//          left = mul(right.data, left);
//      }
//
// Instead, we will define an ordinary function for the one case that
// we've run into in a test effect so far:
//
void mulAssign(inout vec2 left, mat2 right)
{
    left = mul(right.data, left);
}

float fract(float value)
{
    return frac(value);
}



// a * b in GLSL!
// float2 mul(float2 b, mat2 a) {
//     return mul(b, a.data);
// }

float2 operator*(mat2 a, float2 b)
{
    return mul(b, a.data); // have to be flipped!
}

float2 operator*(float2 a, mat2 b)
{
    return mul(b.data, a);
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
// float3 mul(float3 b, mat3 a) {
//     return mul(b, a.data);
// }

float3 operator*(mat3 a, float3 b)
{
    return mul(b, a.data); // have to be flipped!
}



vec2 dFdx(vec2 v)
{
    return ddx(v);
}
vec2 dFdy(vec2 v)
{
    return ddy(v);
}

float atan(float y, float x)
{
    return atan2(y, x);
}


typealias sampler2D = TextureSamplerPair;



/// ----------------------------------------------------------------
/// User code will go here ------------------------------------------------

