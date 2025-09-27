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
    float4 iMouse;
    float2 iResolution;
    float iTime;
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

float3 operator*(mat3 a, float3 b)
{
    return mul(b, a); // have to be flipped!
}

typealias sampler2D = TextureSamplerPair;



/// ----------------------------------------------------------------
/// User code will go here ------------------------------------------------


// https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang
// ^ Adapted this removing the interface, which was a lot of bloat for not a lot of functionality
// the final slang will be a sandwich with shadertoy1.slang + user code (beginning of shadertoy2) + shadertoy2.slang

struct ShadertoyWrapper {
    // User code will go here and provide the impl for:
    // void mainImage(out float4 fragColor, in float2 fragCoord); 

    // ---------------- User Code Starts Here ----------------
    #define T (iTime*.5)
#define P(z) (vec3((cos((z) * .2) * .3) * 10., \
                   (cos((z) * .15) * .3) * 10., (z)))
#define R(a) mat2(cos(a), -sin(a), sin(a), cos(a))
#define N normalize

bool hitRoot = false;

// @Shane
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n ){
  
    n = max((abs(n) - 0.2)*7., 0.001);
    n /= (n.x + n.y + n.z );  
    
	return (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
}

// @Shane:
// Commutative smooth minimum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smin(float a, float b, float k){

   float f = max(0., 1. - abs(b - a)/k);
   return min(a, b) - k*.25*f*f;
}


float fractal(vec3 p) {
    float s,w,l;

    p.y *= .4;
    p += cos(p.yzx*12.)*.07;
    p.x -= 1.6;
    for (s=0.,w=.7; s++ < 7.; p *= l, w *= l )
        p  = abs(sin(p))-1.,
        l = 1.5/dot(p,p);

    return length(p)/w-.0005;
}

float light(vec3 p) {
    p.xy -= P(p.z).xy;
    float t = iTime*.5;
    return length(p-vec3(
            tanh(cos(t*1.8)*2.) * .25 ,
            tanh(cos(t*1.)*3.) * .3,
            4.+T+tanh(cos(t*.4)*4.)*1.
        ))-.25;
}

float root(vec3 p){
    p.xy *= R(sin(p.z*2.)*.1 + sin(p.z*1.)*.3);
    p = mod(p, 1.) - .5;
    return .6*length(p.xy) - .01;
}

float map(vec3 p) {
    float s, r;
    r = root(p);
    
    p.xy -= P(p.z).xy;
    s = smin(fractal(p), 2. - abs(p.y), .8);
    hitRoot = r < s;
    s = min(s, r);
    return s;
}

// @iq
float AO(in vec3 pos, in vec3 nor) {
	float sca = 1.9, occ = 0.0;
    for( int i=0; i<5; i++ ){
        float hr = 0.01 + float(i)*0.5/4.0;        
        float dd = map(nor * hr + pos);
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp( 1.0 - occ, 0.0, 1.0 );    
}

void mainImage(out vec4 o, in vec2 u) {
    float s=.1,d=0.,i=0.;
    vec3  r = iResolution;
    u = (u-r.xy/2.)/r.y;
    vec3  e = vec3(.001,0,0),
          p = P(T),ro=p,
          Z = N( P(T+1.) - p),
          X = N(vec3(Z.z,0,-Z)),
          D = vec3(R(sin(T*.2)*.4)*u, 1) 
             * mat3(-X, cross(X, Z), Z);

    o = vec4(0);
    for(;i++ < 128.;)
        p = ro + D * d * .65,
        d += s = map(p),
        o += 4e1*vec4(6,9,1,0) + 1./max(s, .001);
        
    r =  N(map(p) - vec3(map(p-e.xyy), 
                         map(p-e.yxy), 
                         map(p-e.yyx)));
    
    o.rgb *= tex3D(iChannel0, p*.5, r);
    o *= o;
    o *= hitRoot ? .2*vec4(4,3,2,0) : vec4(1);
    o *= AO(p, r);
    o = tanh(sqrt(d * o / pow(abs(light(p)), 3.) / 5e9 ));

}
    // ---------------- User Code Ends Here ------------------

    static This getDefault() { return This(); }
}

ShadertoyWrapper g_shadertoy;

uniform RWTexture2D<float4> g_output;

// By defining an interface for image shader effects, we
// have been able to decouple the code for the effects
// themselves from the code for their execution contexts.
// A key benefit of that decoupling is that we can introduce
// both new effects and new execution contexts in a modular
// fashion.
//
// For example, we can easily define a compute shader for
// executing an image shader effect:
//
[shader("compute")]
[numthreads(16, 1, 1)]
void compute_main(
    uint3 sv_dispatchThreadID : SV_DispatchThreadID)
{
    // The operations required to set up and execute
    // the user-defined effect are similar to what
    // they were for the fragment shader.
    //
    ShadertoyWrapper toy = ShadertoyWrapper.getDefault();

    float2 fragCoord = float2(sv_dispatchThreadID.xy);
    fragCoord.y = iResolution.y - fragCoord.y - 1; // Flip Y!!!
    float4 fragColor = 0;
    g_shadertoy.mainImage(fragColor, fragCoord);

    // The main difference is that we now write the
    // output color explicitly to an image pixel
    // instead of relying on the rasterization pipeline.
    //
    g_output[sv_dispatchThreadID.xy] = fragColor;
}

