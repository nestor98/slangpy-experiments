// https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang
// ^ Adapted this removing the interface, which was a lot of bloat for not a lot of functionality
// the final slang will be a sandwich with shadertoy1.slang + user code (beginning of shadertoy2) + shadertoy2.slang


struct ShadertoyWrapper {
    // SPECIFIC DEFINES HERE: !!!!!! WARNING: THEY MAKE THE PROGRAM CRASH RIGHT NOW
    static float tube_end_x = 10.;

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
        // x_range = float2(-10.0, 10.);
        tube_end_x = -10.;

        min_f = 0.9;
        max_f = 2.; // min and max freq
        min_twist = -.34;
        max_twist = -.1;

        wl_idx = float[11](
            0.2, // blue
            0.9, // red
            0.4, // green
            .5,  // yellow
            0.6, // orange
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

uniform RWTexture2D<float4> g_output;

[shader("compute")]
[numthreads(16, 1, 1)]
[mutating]
void compute_main(
    uint3 sv_dispatchThreadID : SV_DispatchThreadID)
{
    ShadertoyWrapper shadertoy = ShadertoyWrapper();
    float2 fragCoord = float2(sv_dispatchThreadID.xy);
    // fragCoord.y = iResolution.y - fragCoord.y - 1; // Flip Y!!!
    float4 fragColor = 0;
    shadertoy.mainImage(fragColor, fragCoord);

    uint2 idx = uint2(sv_dispatchThreadID.x, uint(iResolution.y) - 1 - sv_dispatchThreadID.y);
    g_output[idx] = fragColor;
}
