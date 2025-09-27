// https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang
// ^ Adapted this removing the interface, which was a lot of bloat for not a lot of functionality
// the final slang will be a sandwich with shadertoy1.slang + user code (beginning of shadertoy2) + shadertoy2.slang

struct ShadertoyWrapper {
    // SPECIFIC DEFINES HERE: !!!!!! WARNING: THEY MAKE THE PROGRAM CRASH RIGHT NOW
    static float2 x_range = float2(0, 10.);

    static float wl_idx[11];

    static float min_f = 0.9, max_f = 2.; // min and max freq
    static float min_twist = -.34, max_twist = -.1;


    // User code will go here and provide the impl for:
    // void mainImage(out float4 fragColor, in float2 fragCoord);

    // ---------------- User Code Starts Here ----------------
    {USER_CODE}
    // ---------------- User Code Ends Here ------------------

    __init() {
        // SPECIFIC INITIALIZATION HERE:
        x_range = float2(0, 10.);

        min_f = 0.9;
        max_f = 2.; // min and max freq
        min_twist = -.34;
        max_twist = -.1;

        wl_idx = float[11](
            0.01, // blue
            0.99, // red
            0.33, // green
            .5,   // yellow
            0.6,  // orange
            0.1,
            0.2,
            0.4,
            0.7,
            0.65,
            0.85 // red
        );
    }
    static This getDefault() { return This(); }
}

ShadertoyWrapper g_shadertoy;

uniform RWTexture2D<float4> g_output;

[shader("compute")]
[numthreads(16, 1, 1)]
[mutating]
void compute_main(
    uint3 sv_dispatchThreadID : SV_DispatchThreadID)
{
    float2 fragCoord = float2(sv_dispatchThreadID.xy);
    fragCoord.y = iResolution.y - fragCoord.y - 1; // Flip Y!!!
    float4 fragColor = 0;
    g_shadertoy.mainImage(fragColor, fragCoord);
    
    g_output[sv_dispatchThreadID.xy] = fragColor;
}

