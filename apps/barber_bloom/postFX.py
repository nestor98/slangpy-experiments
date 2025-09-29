# /***************************************************************************
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
#  **************************************************************************/

# This file is an adaptation of the falcor .h file at https://github.com/NVIDIAGameWorks/Falcor/blob/eb540f6748774680ce0039aaf3ac9279266ec521/Source/RenderPasses/SimplePostFX/

import slangpy as spy
import numpy as np

from slangpy.types import call_id


class PostFX:
    """
    A wrapper around the postFX.slang file using a slangpy module. The ui uses spy.ui, which right now is super ugly... could be switched for another one
    Usage: 
    init: postfx = PostFX(device)
          postfx.setup_ui(parent)  # parent is a spy.ui.Widget

          loop:
            self.window.process_events()
            self.ui.process_events()
            

            texture2D_hdr = ... # previous passes

            postfx.set(...)  # override properties as needed
            
            postfx.execute(pRenderContext, texture2D_hdr, output)

    """
    def __init__(self, device):
        self.kNumLevels = 8  # Constant for the number of levels

        # Selected output size
        # self.mOutputSizeSelection = "Default"  # Placeholder for RenderPassHelpers::IOSize
        # Output size in pixels when 'Fixed' size is selected
        self.mFixedOutputSize = spy.uint2(512, 512)

        self.mpDownsamplePass = None
        self.mpUpsamplePass = None
        self.mpPostFXPass = None

        # Image pyramid, fine to coarse, full res down in steps of 4x (16x area)
        self.mpPyramid = [None] * (self.kNumLevels + 1)
        self.mpLinearSampler = None

        # Wipe across to see the effect without fx. 0<=all effect, 1>= disabled
        self.mWipe = 0.0
        # Enable the entire pass
        self.mEnabled = True

        # Amount of bloom
        self.mBloomAmount = 0.01
        # Some extra controls for bloom:
        self.mBloomStrength = 1.0
        self.mBloomThreshold = 1.0
        self.mBloomClamp = 100.0
        self.mBloomSize = 1.0

        # How much of a 6-pointed star to add to the bloom kernel
        self.mStarAmount = 0.0
        # Angle of star rays
        self.mStarAngle = 0.1
        # Amount of circuit vignetting
        self.mVignetteAmount = 0.0
        # Amount of radial chromatic aberration
        self.mChromaticAberrationAmount = 0.0
        # Amount of barrel distortion
        self.mBarrelDistortAmount = 0.0
        # Saturation amount for shadows, midtones, and highlights
        self.mSaturationCurve = spy.float3(1.0, 1.0, 1.0)
        # Color offset, tints shadows
        self.mColorOffset = spy.float3(0.5, 0.5, 0.5)
        # Color scale, tints highlights
        self.mColorScale = spy.float3(0.5, 0.5, 0.5)
        # Color power (gamma), tints midtones
        self.mColorPower = spy.float3(0.5, 0.5, 0.5)

        # The above colors are also offered as scalars for ease of UI and also to set negative colors

        # Luma offset, crushes shadows if negative
        self.mColorOffsetScalar = 0.0
        # Luma scale, effectively another exposure control
        self.mColorScaleScalar = 0.0
        # Luma power, i.e., a gamma curve
        self.mColorPowerScalar = 0.0

        self.mTonemap = True

        self.device = device
        self.module = spy.Module.load_from_file(device, "postFX.mod.slang")

    def set(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Unknown property '{key}' in SimplePostFX properties.")

    def setup_ui(self, parent: spy.ui.Widget):
        widget = spy.ui.Group(parent, "PostFX Settings")
        spy.ui.CheckBox(widget, "Enable post fx", value=self.mEnabled, callback=lambda v: setattr(self, 'mEnabled', v))
        spy.ui.CheckBox(widget, "Tonemap", value=self.mTonemap, callback=lambda v: setattr(self, 'mTonemap', v))
        self.mWipe = spy.ui.SliderFloat(widget, "Wipe", value=self.mWipe, min=0.0, max=1.0)
        lens_fx_group = spy.ui.Group(widget, "Lens FX")

        bloom_group = spy.ui.Group(lens_fx_group, "Bloom")
        self.mBloomAmount = spy.ui.SliderFloat(bloom_group, "Bloom", value=self.mBloomAmount, min=0.0, max=1.0, flags=spy.ui.SliderFlags.logarithmic)
        self.mBloomStrength = spy.ui.SliderFloat(bloom_group, "Bloom Strength", value=self.mBloomStrength, min=0.0, max=10.0, flags=spy.ui.SliderFlags.logarithmic)
        self.mBloomThreshold = spy.ui.SliderFloat(bloom_group, "Bloom Threshold", value=self.mBloomThreshold, min=0.0, max=1000.0, flags=spy.ui.SliderFlags.logarithmic)
        self.mBloomClamp = spy.ui.SliderFloat(bloom_group, "Bloom Clamp", value=self.mBloomClamp, min=1.0, max=1000.0, flags=spy.ui.SliderFlags.logarithmic)
        self.mBloomSize = spy.ui.SliderFloat(bloom_group, "Bloom Size", value=self.mBloomSize, min=0.0, max=5.0, flags=spy.ui.SliderFlags.logarithmic)

        self.mStarAmount = spy.ui.SliderFloat(lens_fx_group, "Bloom Star", value=self.mStarAmount, min=0.0, max=1.0, flags=spy.ui.SliderFlags.logarithmic)
        self.mStarAngle = spy.ui.SliderFloat(lens_fx_group, "Star Angle", value=self.mStarAngle, min=0.0, max=1.0)
        self.mVignetteAmount = spy.ui.SliderFloat(lens_fx_group, "Vignette", value=self.mVignetteAmount, min=0.0, max=1.0)
        self.mChromaticAberrationAmount = spy.ui.SliderFloat(lens_fx_group, "Chromatic Aberration", value=self.mChromaticAberrationAmount, min=0.0, max=1.0)
        self.mBarrelDistortAmount = spy.ui.SliderFloat(lens_fx_group, "Barrel Distortion", value=self.mBarrelDistortAmount, min=0.0, max=1.0)
        def reset_lens_fx():
            self.mBloomAmount.value = 0.01
            self.mBloomStrength.value = 1.0
            self.mBloomThreshold.value = 1.0
            self.mBloomSize.value = 1.0
            self.mBloomClamp.value = 100.0
            self.mStarAmount.value = 0.0
            self.mStarAngle.value = 0.1
            self.mVignetteAmount.value = 0.0
            self.mChromaticAberrationAmount.value = 0.0
            self.mBarrelDistortAmount.value = 0.0
        spy.ui.Button(lens_fx_group, "reset this group", callback=reset_lens_fx)

        saturation_group = spy.ui.Group(widget, "Saturation")
        self.mSaturationCurve = spy.ui.SliderFloat3(saturation_group, "Saturation Curve", value=self.mSaturationCurve, min=0.0, max=2.0)
        def reset_saturation():
            self.mSaturationCurve.value = spy.float3(1.0, 1.0, 1.0)
        spy.ui.Button(saturation_group, "reset this group", callback=reset_saturation)

        luma_group = spy.ui.Group(widget, "Offset/Power/Scale (luma)")
        self.mColorOffsetScalar = spy.ui.SliderFloat(luma_group, "Luma Offset (Shadows)", value=self.mColorOffsetScalar, min=-1.0, max=1.0)
        self.mColorPowerScalar = spy.ui.SliderFloat(luma_group, "Luma Power (Midtones)", value=self.mColorPowerScalar, min=-1.0, max=1.0)
        self.mColorScaleScalar = spy.ui.SliderFloat(luma_group, "Luma Scale (Hilights)", value=self.mColorScaleScalar, min=-1.0, max=1.0)
        def reset_luma():
            self.mColorOffsetScalar.value = 0.
            self.mColorPowerScalar.value = 0.0
            self.mColorScaleScalar.value = 0.0

        spy.ui.Button(luma_group, "reset this group", callback=reset_luma)

        self.capture_renderdoc_frame = False
        spy.ui.Button(widget, "Capture to RenderDoc", callback=lambda: setattr(self, 'capture_renderdoc_frame', True))

    def get_bloom_uniforms_dict(self, resolution=spy.uint2(512, 512)):
        return {
            "gResolution"  : resolution,
            "gInvRes": spy.float2(1.0 / float(resolution.x), 1.0 / float(resolution.y)),
            "gBloomAmount": self.mBloomAmount.value,
            "gBloomStrength": self.mBloomStrength.value,
            "gBloomThreshold": self.mBloomThreshold.value,
            "gBloomClamp": self.mBloomClamp.value,
            # "gStarAmount": self.mStarAmount.value,
            # "gStarAngle": self.mStarAngle.value,
        }

    def get_uniforms_dict(self, resolution=spy.uint2(512, 512)):
        satcurve = spy.float3(self.mSaturationCurve.value)
        # Fit a quadratic through the 3 points
        satcurve.y -= satcurve.x
        satcurve.z -= satcurve.x
        A = 2.0 * satcurve.z - 4.0 * satcurve.y
        B = satcurve.z - A
        C = satcurve.x
        gSaturationCurve = spy.float3(A, B, C)

        # gBloomed = self.mpPyramid[0] if self.mBloomAmount.value > 0.0 else None
        
        barrel =  0.125 * self.mBarrelDistortAmount.value

        return {
            "gWipe": self.mWipe.value * resolution.x,
            "gResolution"  : resolution,
            "gInvRes": spy.float2(1.0 / float(resolution.x), 1.0 / float(resolution.y)),
            # "gEnabled": self.mEnabled,
            # "gBloomAmount": self.mBloomAmount.value,
            "gSaturationCurve": gSaturationCurve,
            # "gStarAmount": self.mStarAmount.value,
            # "gStarAngle": self.mStarAngle.value,
            "gVignetteAmount": self.mVignetteAmount.value,
            "gChromaticAberrationAmount": self.mChromaticAberrationAmount.value / 64.,
            # "gBarrelDistortAmount": self.mBarrelDistortAmount.value,
            "gBarrelDistort" : spy.float2(1. / (1. + 4. * barrel), barrel),
            "gSaturationCurve": gSaturationCurve,
            "gColorOffset": self.mColorOffset + self.mColorOffsetScalar.value - 0.5,
            "gColorScale": self.mColorScale * 2**(1.0 + 2.0 * self.mColorScaleScalar.value),
            "gColorPower": spy.float3(2**(3.0 * (0.5 - self.mColorPower.x - self.mColorPowerScalar.value))),        
            "gTonemap": self.mTonemap
        }
    

    def prepare(self, width, height):
        for res in range(self.kNumLevels + 1):
            if self.mBloomAmount.value <= 0.0:
                self.mpPyramid[res] = None
            else:
                w = max(1, width >> res)
                h = max(1, height >> res)
                if self.mpPyramid[res] is None or self.mpPyramid[res].width != w or self.mpPyramid[res].height != h:
                    self.mpPyramid[res] = self.device.create_texture(
                        width=w,
                        height=h,
                        format=spy.Format.rgba16_float,
                        mip_count=1,
                        # None,
                        usage=spy.TextureUsage.shader_resource | spy.TextureUsage.unordered_access
                    )

    def execute(self, command_encoder, pSrc, pDst):
        """
        command_encoder: command encoder
        pSrc: slangpy.Texture2D
        pDst: slangpy.Texture2D
        """
        assert pSrc is not None and pDst is not None

        if self.capture_renderdoc_frame:
            ok = spy.renderdoc.start_frame_capture(device=self.device)
            if ok:
                spy.Logger().info("Capturing to RenderDoc")
            else:
                spy.Logger().error("Failed to start RenderDoc capture (make sure you launched this app from RenderDoc!)")
                self.capture_renderdoc_frame = False

        # Issue error and disable pass if I/O size doesn't match. The user can hit continue and fix the config or abort.
        if self.mEnabled and (pSrc.width != pDst.width or pSrc.height != pDst.height):
            spy.Logger().error("SimplePostFX I/O sizes don't match. The pass will be disabled.")
            self.mEnabled = False
        resolution = spy.uint2(pSrc.width, pSrc.height)

        # if we have 'identity' settings, we can just copy input to output
        if not self.mEnabled or self.mWipe.value >= 1.0 or (
            self.mBloomAmount.value <= 0.0 and
            self.mVignetteAmount.value <= 0.0 and
            self.mChromaticAberrationAmount.value <= 0.0 and
            self.mBarrelDistortAmount.value <= 0.0 and
            all(self.mSaturationCurve.value == spy.float3(1.0)) and
            all(self.mColorOffset == spy.float3(0.5)) and
            all(self.mColorScale == spy.float3(0.5)) and
            all(self.mColorPower == spy.float3(0.5)) and
            self.mColorOffsetScalar.value <= 0.0 and
            self.mColorScaleScalar.value <= 0.0 and
            self.mColorPowerScalar.value <= 0.0
        ):
            # wipe is all the way across, which corresponds to no effect
            # pRenderContext.blit(pDst, pSrc)
            return False

        self.prepare(resolution.x, resolution.y)

        sampler = self.device.create_sampler(
            address_u=spy.TextureAddressingMode.clamp_to_border, address_v=spy.TextureAddressingMode.clamp_to_border
            , address_w=spy.TextureAddressingMode.clamp_to_border,
            border_color=spy.float4(0.0, 0.0, 0.0, 0.0)
        ) 
        # downsample:
        if self.mBloomAmount.value > 0.0:
            assert self.mpPyramid[0] is not None
            for level in range(self.kNumLevels):
                assert self.mpPyramid[level + 1] is not None
                res = spy.uint2(max(1, resolution.x >> (level + 1)), max(1, resolution.y >> (level + 1)))
                invres = spy.float2(1.0 / res.x, 1.0 / res.y)
                self.module.downsample(pixel=spy.call_id((res.y, res.x)),# spy.grid(res.x, res.y),  # spy.int2(10,10),#
                                       _append_to=command_encoder,
                                       cb={"gResolution": res,
                                           "gInvRes": invres,
                                           "gBloomThreshold": self.mBloomThreshold.value,
                                           "gBloomRadius": self.mBloomSize.value
                                           },
                                       gSampler=sampler,
                                       gSrc=self.mpPyramid[level] if level>0 else pSrc,
                                       gDst=self.mpPyramid[level + 1]
                                       )
                # command_encoder.global_barrier()
            # upsample:
            for level in reversed(range(self.kNumLevels)):
                res = spy.uint2(max(1, resolution.x >> level), max(1, resolution.y >> level))
                invres = spy.float2(1.0 / res.x, 1.0 / res.y)
                wantStar = level == 1 or level == 2
                ang = self.mStarAngle.value
                starDir1 = spy.float2(np.sin(ang), np.cos(ang)) * invres * 2.0 if wantStar else spy.float2(0.0, 0.0)
                ang += float(np.pi) / 3.0
                starDir2 = spy.float2(np.sin(ang), np.cos(ang)) * invres * 2.0 if wantStar else spy.float2(0.0, 0.0)
                ang += float(np.pi) / 3.0
                starDir3 = spy.float2(np.sin(ang), np.cos(ang)) * invres * 2.0 if wantStar else spy.float2(0.0, 0.0)
                self.module.upsample(pixel=spy.call_id((res.y, res.x)),# spy.grid(res.x, res.y),  # spy.int2(10,10),#
                                     _append_to=command_encoder,
                                     cb={"gResolution": res, "gInvRes": invres,
                                         "gBloomAmount": self.mBloomAmount.value,
                                         "gBloomStrength": self.mBloomStrength.value,
                                         "gBloomThreshold": self.mBloomThreshold.value,
                                         "gBloomClamp": self.mBloomClamp.value,
                                         "gBloomRadius": self.mBloomSize.value,
                                         "gStar": self.mStarAmount.value if wantStar else 0.0,
                                         "gStarDir1": starDir1, "gStarDir2": starDir2, "gStarDir3": starDir3, "gInPlace": level > 0},
                                     gSampler=sampler,
                                     gSrc=pSrc,
                                     gBloomed=self.mpPyramid[level + 1],
                                     gDst=self.mpPyramid[level] #if level > 0 else pSrc
                                     )
                # command_encoder.global_barrier()
                
        # -----------------------------------------------
        # upsample
        # print("res: ", resolution)
        self.module.runPostFX(pixel=spy.call_id((resolution.y, resolution.x)),# spy.grid(resolution.x, resolution.y),  # spy.int2(10,10),#
                              _append_to=command_encoder,
                              cb=self.get_uniforms_dict(resolution),
                              gSrc=pSrc,
                              gBloomed=self.mpPyramid[0] if self.mBloomAmount.value > 0.0 else pSrc,
                              gSampler=sampler,
                              gDst=pDst
                              )
        
        if self.capture_renderdoc_frame:
            spy.renderdoc.end_frame_capture()
            self.capture_renderdoc_frame = False
            spy.Logger().info("Captured frame. Go to [File>Attach to Running Instance] in RenderDoc to see it.")

        return True 
    