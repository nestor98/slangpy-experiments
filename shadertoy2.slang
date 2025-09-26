// https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang
// ^ Adapted this removing the interface, which was a lot of bloat for not a lot of functionality
// the final slang will be a sandwich with shadertoy1.slang + user code (beginning of shadertoy2) + shadertoy2.slang

struct ShadertoyWrapper {
    // User code will go here and provide the impl for:
    // void mainImage(out float4 fragColor, in float2 fragCoord); 

    // ---------------- User Code Starts Here ----------------
    {USER_CODE}
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
[mutating]
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

