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
        self.mBloomAmount = 0.0
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
        self.mEnabled = spy.ui.CheckBox(widget, "Enable post fx", value=self.mEnabled, callback=lambda v: setattr(self, 'mEnabled', v))
        self.mWipe = spy.ui.SliderFloat(widget, "Wipe", value=self.mWipe, min=0.0, max=1.0)
        lens_fx_group = spy.ui.Group(widget, "Lens FX")
        self.mBloomAmount = spy.ui.SliderFloat(lens_fx_group, "Bloom", value=self.mBloomAmount, min=0.0, max=1.0)
        self.mStarAmount = spy.ui.SliderFloat(lens_fx_group, "Bloom Star", value=self.mStarAmount, min=0.0, max=1.0)
        self.mStarAngle = spy.ui.SliderFloat(lens_fx_group, "Star Angle", value=self.mStarAngle, min=0.0, max=1.0, flags=spy.ui.SliderFlags.logarithmic)
        self.mVignetteAmount = spy.ui.SliderFloat(lens_fx_group, "Vignette", value=self.mVignetteAmount, min=0.0, max=1.0)
        self.mChromaticAberrationAmount = spy.ui.SliderFloat(lens_fx_group, "Chromatic Aberration", value=self.mChromaticAberrationAmount, min=0.0, max=1.0)
        self.mBarrelDistortAmount = spy.ui.SliderFloat(lens_fx_group, "Barrel Distortion", value=self.mBarrelDistortAmount, min=0.0, max=1.0)
        def reset_lens_fx():
            self.mBloomAmount.value = 0.0
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


    def get_uniforms_dict(self):
        return {
            "gWipe": self.mWipe.value,
            # "gEnabled": self.mEnabled,
            "gBloomAmount": self.mBloomAmount.value,
            # "gStarAmount": self.mStarAmount.value,
            # "gStarAngle": self.mStarAngle.value,
            "gVignetteAmount": self.mVignetteAmount.value,
            "gChromaticAberrationAmount": self.mChromaticAberrationAmount.value,
            # "gBarrelDistortAmount": self.mBarrelDistortAmount.value,
            "gSaturationCurve": self.mSaturationCurve.value,
            "gColorOffset": self.mColorOffset + spy.float3(self.mColorOffsetScalar.value) - 0.5,
            "gColorScale": self.mColorScale * 2**(1.0 + 2.0 * self.mColorScaleScalar.value),
            "gColorPower": spy.float3(2**(3.0 * (0.5 - self.mColorPower.x - self.mColorPowerScalar.value))),
        }

        # i dont think there are colorpickers in python :(
        # color_group = spy.ui.Group(widget, "Offset/Power/Scale (color)", True)
        # if spy.ui.Button(color_group, "reset##1"):
        #     self.mColorOffset = spy.float3(0.5, 0.5, 0.5)
        # spy.ui.ColorPicker(color_group, "Color Offset (Shadows)", self.mColorOffset, True)

        # if spy.ui.Button(color_group, "reset##2"):
        #     self.mColorPower = spy.float3(0.5, 0.5, 0.5)
        # spy.ui.ColorPicker(color_group, "Color Power (Midtones)", self.mColorPower, True)

        # if spy.ui.Button(color_group, "reset##3"):
        #     self.mColorScale = spy.float3(0.5, 0.5, 0.5)
        # spy.ui.ColorPicker(color_group, "Color Scale (Hilights)", self.mColorScale, True)

        # if spy.ui.Button(color_group, "reset this group"):
        #     self.mColorOffset = spy.float3(0.5, 0.5, 0.5)
        #     self.mColorPower = spy.float3(0.5, 0.5, 0.5)
        #     self.mColorScale = spy.float3(0.5, 0.5, 0.5)

    def prepare(self, pRenderContext, width, height):
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
                        # 1,
                        # 1,
                        # None,
                        usage=spy.TextureUsage.shader_resource | spy.TextureUsage.unordered_access
                    )

    def execute(self, pRenderContext, pSrc, pDst):
        """
        pRenderContext: slangpy.RenderContext
        pSrc: slangpy.Texture2D
        pDst: slangpy.Texture2D
        """
        assert pSrc is not None and pDst is not None

        # Issue error and disable pass if I/O size doesn't match. The user can hit continue and fix the config or abort.
        if self.mEnabled and (pSrc.width != pDst.width or pSrc.height != pDst.height):
            spy.Logger().error("SimplePostFX I/O sizes don't match. The pass will be disabled.")
            self.mEnabled = False
        resolution = spy.uint2(pSrc.width, pSrc.height)

        # if we have 'identity' settings, we can just copy input to output
        # if not self.mEnabled or self.mWipe.value >= 1.0 or (
        #     self.mBloomAmount.value <= 0.0 and
        #     self.mChromaticAberrationAmount.value <= 0.0 and
        #     self.mBarrelDistortAmount.value <= 0.0 and
        #     all(self.mSaturationCurve.value == spy.float3(1.0)) and
        #     all(self.mColorOffset == spy.float3(0.5)) and
        #     all(self.mColorScale == spy.float3(0.5)) and
        #     all(self.mColorPower == spy.float3(0.5)) and
        #     self.mColorOffsetScalar.value <= 0.0 and
        #     self.mColorScaleScalar.value <= 0.0 and
        #     self.mColorPowerScalar.value <= 0.0
        # ):
        #     # wipe is all the way across, which corresponds to no effect
        #     # pRenderContext.blit(pDst, pSrc)
        #     return False

        self.prepare(pRenderContext, resolution.x, resolution.y)
        # if self.mBloomAmount > 0.0:
        #     # Downsampling
        #     var = self.mpDownsamplePass.get_root_var()
        #     var["gLinearSampler"] = self.module.get_device().create_sampler(spy.Sampler.Desc().set_filter_mode(
        #         spy.TextureFilteringMode.Linear, spy.TextureFilteringMode.Linear, spy.TextureFilteringMode.Point).set_addressing_mode(
        #         spy.TextureAddressingMode.Border, spy.TextureAddressingMode.Border, spy.TextureAddressingMode.Border))
        #     for level in range(self.kNumLevels):
        #         res = spy.uint2(max(1, resolution.x >> (level + 1)), max(1, resolution.y >> (level + 1)))
        #         invres = spy.float2(1.0 / res.x, 1.0 / res.y)
        #         var["PerFrameCB"]["gResolution"] = res
        #         var["PerFrameCB"]["gInvRes"] = invres
        #         var["gSrc"] = self.mpPyramid[level] if level else pSrc
        #         var["gDst"] = self.mpPyramid[level + 1]
        #         self.mpDownsamplePass.execute(pRenderContext, spy.uint3(res.x, res.y, 1))
        #     # Upsampling
        #     var = self.mpUpsamplePass.get_root_var()
        #     var["gLinearSampler"] = self.module.get_device().create_sampler(spy.Sampler.Desc().set_filter_mode(
        #         spy.TextureFilteringMode.Linear, spy.TextureFilteringMode.Linear, spy.TextureFilteringMode.Point).set_addressing_mode(
        #         spy.TextureAddressingMode.Border, spy.TextureAddressingMode.Border, spy.TextureAddressingMode.Border))
        #     var["PerFrameCB"]["gBloomAmount"] = self.mBloomAmount
        #     var["gSrc"] = pSrc
        #     for level in reversed(range(self.kNumLevels)):
        #         res = spy.uint2(max(1, resolution.x >> level), max(1, resolution.y >> level))
        #         invres = spy.float2(1.0 / res.x, 1.0 / res.y)
        #         var["PerFrameCB"]["gResolution"] = res
        #         var["PerFrameCB"]["gInvRes"] = invres
        #         wantStar = level == 1 or level == 2
        #         var["PerFrameCB"]["gStar"] = self.mStarAmount if wantStar else 0.0
        #         if wantStar:
        #             ang = self.mStarAngle
        #             var["PerFrameCB"]["gStarDir1"] = spy.float2(spy.sin(ang), spy.cos(ang)) * invres * 2.0
        #             ang += float(spy.pi) / 3.0
        #             var["PerFrameCB"]["gStarDir2"] = spy.float2(spy.sin(ang), spy.cos(ang)) * invres * 2.0
        #             ang += float(spy.pi) / 3.0
        #             var["PerFrameCB"]["gStarDir3"] = spy.float2(spy.sin(ang), spy.cos(ang)) * invres * 2.0
        #         var["gBloomed"] = self.mpPyramid[level + 1]
        #         var["gDst"] = self.mpPyramid[level]
        #         # for most levels, we update the pyramid in place. for the last step, we read
        #         # from the original source since we did not compute it in the downsample passes.
        #         var["PerFrameCB"]["gInPlace"] = level > 0
        #         self.mpUpsamplePass.execute(pRenderContext, spy.uint3(res.x, res.y, 1))
        # -----------------------------------------------
        # upsample
        print(call_id())
        self.module.runPostFX(pixel=spy.int2(call_id()),# spy.grid(resolution.x, resolution.y),  # spy.int2(10,10),#
                            #   _append_to=pRenderContext,
                              cb=self.get_uniforms_dict(),
                              gSrc=pSrc,
                              gBloomed=self.mpPyramid[0] if self.mBloomAmount.value > 0.0 else pSrc,
                              gSampler=self.device.create_sampler(),
                              gDst=pDst
                                #   {"PerFrameCB": self.get_uniforms_dict(),
                                #     # {"gLinearSampler": self.device.create_sampler(spy.Sampler.Desc().set_filter_mode(
                                #     #     spy.TextureFilteringMode.Linear, spy.TextureFilteringMode.Linear, spy.TextureFilteringMode.Point).set_addressing_mode(
                                #     #     spy.TextureAddressingMode.Border, spy.TextureAddressingMode.Border, spy.TextureAddressingMode.Border))},
                                #     "gSrc": pSrc,
                                #     "gDst": pDst,
                                #     "gBloomed": self.mpPyramid[0] if self.mBloomAmount.value > 0.0 else pSrc}
                              

                              )
        
        return True 
# CPP: 
# /***************************************************************************
#  # Copyright (c) 2015-23, NVIDIA CORPORATION. All rights reserved.
#  #
#  # Redistribution and use in source and binary forms, with or without
#  # modification, are permitted provided that the following conditions
#  # are met:
#  #  * Redistributions of source code must retain the above copyright
#  #    notice, this list of conditions and the following disclaimer.
#  #  * Redistributions in binary form must reproduce the above copyright
#  #    notice, this list of conditions and the following disclaimer in the
#  #    documentation and/or other materials provided with the distribution.
#  #  * Neither the name of NVIDIA CORPORATION nor the names of its
#  #    contributors may be used to endorse or promote products derived
#  #    from this software without specific prior written permission.
#  #
#  # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
#  # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  # PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
#  # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#  # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#  # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
#  # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  **************************************************************************/
# #include "SimplePostFX.h"

# namespace
# {
# const char kSrc[] = "src";
# const char kDst[] = "dst";

# // Scripting options.
# const char kEnabled[] = "enabled";
# const char kOutputSize[] = "outputSize";
# const char kFixedOutputSize[] = "fixedOutputSize";
# const char kWipe[] = "wipe";
# const char kBloomAmount[] = "bloomAmount";
# const char kStarAmount[] = "starAmount";
# const char kStarAngle[] = "starAngle";
# const char kVignetteAmount[] = "vignetteAmount";
# const char kChromaticAberrationAmount[] = "chromaticAberrationAmount";
# const char kBarrelDistortAmount[] = "barrelDistortAmount";
# const char kSaturationCurve[] = "saturationCurve";
# const char kColorOffset[] = "colorOffset";
# const char kColorScale[] = "colorScale";
# const char kColorPower[] = "colorPower";
# const char kColorOffsetScalar[] = "colorOffsetScalar";
# const char kColorScaleScalar[] = "colorScaleScalar";
# const char kColorPowerScalar[] = "colorPowerScalar";

# const char kShaderFile[] = "RenderPasses/SimplePostFX/SimplePostFX.cs.slang";
# } // namespace

# static void regSimplePostFX(pybind11::module& m)
# {
#     pybind11::class_<SimplePostFX, RenderPass, ref<SimplePostFX>> pass(m, "SimplePostFX");
#     pass.def_property(kEnabled, &SimplePostFX::getEnabled, &SimplePostFX::setEnabled);
#     pass.def_property(kWipe, &SimplePostFX::getWipe, &SimplePostFX::setWipe);
#     pass.def_property(kBloomAmount, &SimplePostFX::getBloomAmount, &SimplePostFX::setBloomAmount);
#     pass.def_property(kStarAmount, &SimplePostFX::getStarAmount, &SimplePostFX::setStarAmount);
#     pass.def_property(kStarAngle, &SimplePostFX::getStarAngle, &SimplePostFX::setStarAngle);
#     pass.def_property(kVignetteAmount, &SimplePostFX::getVignetteAmount, &SimplePostFX::setVignetteAmount);
#     pass.def_property(kChromaticAberrationAmount, &SimplePostFX::getChromaticAberrationAmount, &SimplePostFX::setChromaticAberrationAmount);
#     pass.def_property(kBarrelDistortAmount, &SimplePostFX::getBarrelDistortAmount, &SimplePostFX::setBarrelDistortAmount);
#     pass.def_property(kSaturationCurve, &SimplePostFX::getSaturationCurve, &SimplePostFX::setSaturationCurve);
#     pass.def_property(kColorOffset, &SimplePostFX::getColorOffset, &SimplePostFX::setColorOffset);
#     pass.def_property(kColorScale, &SimplePostFX::getColorScale, &SimplePostFX::setColorScale);
#     pass.def_property(kColorPower, &SimplePostFX::getColorPower, &SimplePostFX::setColorPower);
#     pass.def_property(kColorOffsetScalar, &SimplePostFX::getColorOffsetScalar, &SimplePostFX::setColorOffsetScalar);
#     pass.def_property(kColorScaleScalar, &SimplePostFX::getColorScaleScalar, &SimplePostFX::setColorScaleScalar);
#     pass.def_property(kColorPowerScalar, &SimplePostFX::getColorPowerScalar, &SimplePostFX::setColorPowerScalar);
# }

# extern "C" FALCOR_API_EXPORT void registerPlugin(Falcor::PluginRegistry& registry)
# {
#     registry.registerClass<RenderPass, SimplePostFX>();
#     ScriptBindings::registerBinding(regSimplePostFX);
# }

# SimplePostFX::SimplePostFX(ref<Device> pDevice, const Properties& props) : RenderPass(pDevice)
# {
#     // Deserialize pass from dictionary.
#     for (const auto& [key, value] : props)
#     {
#         if (key == kEnabled)
#             setEnabled(value);
#         else if (key == kOutputSize)
#             mOutputSizeSelection = value;
#         else if (key == kFixedOutputSize)
#             mFixedOutputSize = value;
#         else if (key == kWipe)
#             setWipe(value);
#         else if (key == kBloomAmount)
#             setBloomAmount(value);
#         else if (key == kStarAmount)
#             setStarAmount(value);
#         else if (key == kStarAngle)
#             setStarAngle(value);
#         else if (key == kVignetteAmount)
#             setVignetteAmount(value);
#         else if (key == kChromaticAberrationAmount)
#             setChromaticAberrationAmount(value);
#         else if (key == kBarrelDistortAmount)
#             setBarrelDistortAmount(value);
#         else if (key == kSaturationCurve)
#             setSaturationCurve(value);
#         else if (key == kColorOffset)
#             setColorOffset(value);
#         else if (key == kColorScale)
#             setColorScale(value);
#         else if (key == kColorPower)
#             setColorPower(value);
#         else if (key == kColorOffsetScalar)
#             setColorOffsetScalar(value);
#         else if (key == kColorScaleScalar)
#             setColorScaleScalar(value);
#         else if (key == kColorPowerScalar)
#             setColorPowerScalar(value);
#         else
#             logWarning("Unknown property '{}' in SimplePostFX properties.", key);
#     }

#     Sampler::Desc samplerDesc;
#     samplerDesc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Point);
#     samplerDesc.setAddressingMode(TextureAddressingMode::Border, TextureAddressingMode::Border, TextureAddressingMode::Border);
#     mpLinearSampler = mpDevice->createSampler(samplerDesc);

#     DefineList defines;
#     mpDownsamplePass = ComputePass::create(mpDevice, kShaderFile, "downsample", defines);
#     mpUpsamplePass = ComputePass::create(mpDevice, kShaderFile, "upsample", defines);
#     mpPostFXPass = ComputePass::create(mpDevice, kShaderFile, "runPostFX", defines);
# }

# Properties SimplePostFX::getProperties() const
# {
#     Properties props;
#     props[kEnabled] = getEnabled();
#     props[kOutputSize] = mOutputSizeSelection;
#     if (mOutputSizeSelection == RenderPassHelpers::IOSize::Fixed)
#         props[kFixedOutputSize] = mFixedOutputSize;
#     props[kWipe] = getWipe();
#     props[kBloomAmount] = getBloomAmount();
#     props[kStarAmount] = getStarAmount();
#     props[kStarAngle] = getStarAngle();
#     props[kVignetteAmount] = getVignetteAmount();
#     props[kChromaticAberrationAmount] = getChromaticAberrationAmount();
#     props[kBarrelDistortAmount] = getBarrelDistortAmount();
#     props[kSaturationCurve] = getSaturationCurve();
#     props[kColorOffset] = getColorOffset();
#     props[kColorScale] = getColorScale();
#     props[kColorPower] = getColorPower();
#     props[kColorOffsetScalar] = getColorOffsetScalar();
#     props[kColorScaleScalar] = getColorScaleScalar();
#     props[kColorPowerScalar] = getColorPowerScalar();
#     return props;
# }

# RenderPassReflection SimplePostFX::reflect(const CompileData& compileData)
# {
#     RenderPassReflection reflector;
#     const uint2 sz = RenderPassHelpers::calculateIOSize(mOutputSizeSelection, mFixedOutputSize, compileData.defaultTexDims);

#     reflector.addInput(kSrc, "Source texture").bindFlags(ResourceBindFlags::ShaderResource);
#     ;
#     reflector.addOutput(kDst, "post-effected output texture")
#         .bindFlags(ResourceBindFlags::RenderTarget | ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess)
#         .format(ResourceFormat::RGBA32Float)
#         .texture2D(sz.x, sz.y);
#     return reflector;
# }

# void SimplePostFX::execute(RenderContext* pRenderContext, const RenderData& renderData)
# {
#     auto pSrc = renderData.getTexture(kSrc);
#     auto pDst = renderData.getTexture(kDst);
#     FALCOR_ASSERT(pSrc && pDst);

#     // Issue error and disable pass if I/O size doesn't match. The user can hit continue and fix the config or abort.
#     if (getEnabled() && (pSrc->getWidth() != pDst->getWidth() || pSrc->getHeight() != pDst->getHeight()))
#     {
#         logError("SimplePostFX I/O sizes don't match. The pass will be disabled.");
#         mEnabled = false;
#     }
#     const uint2 resolution = uint2(pSrc->getWidth(), pSrc->getHeight());

#     // if we have 'identity' settings, we can just copy input to output
#     // clang-format off
#     if (getEnabled() == false || getWipe() >= 1.f || (
#         getBloomAmount() == 0.f &&
#         getChromaticAberrationAmount() == 0.f &&
#         getBarrelDistortAmount() == 0.f &&
#         all(getSaturationCurve() == float3(1.f)) &&
#         all(getColorOffset() == float3(0.5f)) &&
#         all(getColorScale() == float3(0.5f)) &&
#         all(getColorPower() == float3(0.5f)) &&
#         getColorOffsetScalar() == 0.f &&
#         getColorScaleScalar() == 0.f &&
#         getColorPowerScalar() == 0.f
#         ))
#     {
#         // wipe is all the way across, which corresponds to no effect
#         pRenderContext->blit(pSrc->getSRV(), pDst->getRTV());
#         return;
#     }
#     // clang-format on

#     preparePostFX(pRenderContext, resolution.x, resolution.y);
#     if (getBloomAmount() > 0.f)
#     {
#         // Downsampling
#         {
#             auto var = mpDownsamplePass->getRootVar();
#             var["gLinearSampler"] = mpLinearSampler;
#             for (int level = 0; level < kNumLevels; ++level)
#             {
#                 uint2 res = {std::max(1u, resolution.x >> (level + 1)), std::max(1u, resolution.y >> (level + 1))};
#                 float2 invres = float2(1.f / res.x, 1.f / res.y);
#                 var["PerFrameCB"]["gResolution"] = res;
#                 var["PerFrameCB"]["gInvRes"] = invres;
#                 var["gSrc"] = level ? mpPyramid[level] : pSrc;
#                 var["gDst"] = mpPyramid[level + 1];
#                 mpDownsamplePass->execute(pRenderContext, uint3(res, 1));
#             }
#         }

#         // Upsampling
#         {
#             auto var = mpUpsamplePass->getRootVar();
#             var["gLinearSampler"] = mpLinearSampler;
#             var["PerFrameCB"]["gBloomAmount"] = getBloomAmount();
#             var["gSrc"] = pSrc;
#             for (int level = kNumLevels - 1; level >= 0; --level)
#             {
#                 uint2 res = {std::max(1u, resolution.x >> level), std::max(1u, resolution.y >> level)};
#                 float2 invres = float2(1.f / res.x, 1.f / res.y);
#                 var["PerFrameCB"]["gResolution"] = res;
#                 var["PerFrameCB"]["gInvRes"] = invres;
#                 bool wantStar = level == 1 || level == 2;
#                 var["PerFrameCB"]["gStar"] = (wantStar) ? getStarAmount() : 0.f;
#                 if (wantStar)
#                 {
#                     float ang = getStarAngle();
#                     var["PerFrameCB"]["gStarDir1"] = float2(std::sin(ang), std::cos(ang)) * invres * 2.f;
#                     ang += float(M_PI) / 3.f;
#                     var["PerFrameCB"]["gStarDir2"] = float2(std::sin(ang), std::cos(ang)) * invres * 2.f;
#                     ang += float(M_PI) / 3.f;
#                     var["PerFrameCB"]["gStarDir3"] = float2(std::sin(ang), std::cos(ang)) * invres * 2.f;
#                 }
#                 var["gBloomed"] = mpPyramid[level + 1];
#                 var["gDst"] = mpPyramid[level];
#                 // for most levels, we update the pyramid in place. for the last step, we read
#                 // from the original source since we did not compute it in the downsample passes.
#                 var["PerFrameCB"]["gInPlace"] = level > 0;
#                 mpUpsamplePass->execute(pRenderContext, uint3(res, 1));
#             }
#         }
#     }

#     {
#         auto var = mpPostFXPass->getRootVar();
#         var["PerFrameCB"]["gResolution"] = resolution;
#         var["PerFrameCB"]["gInvRes"] = float2(1.f / resolution.x, 1.f / resolution.y);
#         var["PerFrameCB"]["gVignetteAmount"] = getVignetteAmount();
#         var["PerFrameCB"]["gChromaticAberrationAmount"] = getChromaticAberrationAmount() * (1.f / 64.f);
#         float barrel = getBarrelDistortAmount() * 0.125f;
#         // scale factor chosen to keep the corners of a 16:9 viewport fixed
#         var["PerFrameCB"]["gBarrelDistort"] = float2(1.f / (1.f + 4.f * barrel), barrel);
#         float3 satcurve = getSaturationCurve();
#         // fit a quadratic thru the 3 points
#         satcurve.y -= satcurve.x;
#         satcurve.z -= satcurve.x;
#         float A = 2.f * satcurve.z - 4.f * satcurve.y;
#         float B = satcurve.z - A;
#         float C = satcurve.x;
#         var["PerFrameCB"]["gSaturationCurve"] = float3(A, B, C);
#         var["PerFrameCB"]["gColorOffset"] = getColorOffset() + getColorOffsetScalar() - 0.5f;
#         var["PerFrameCB"]["gColorScale"] = getColorScale() * std::exp2(1.f + 2.f * getColorScaleScalar());
#         var["PerFrameCB"]["gColorPower"] = exp2(3.f * (0.5f - getColorPower() - getColorPowerScalar()));
#         var["PerFrameCB"]["gWipe"] = mWipe * resolution.x;
#         var["gBloomed"] = getBloomAmount() > 0.f ? mpPyramid[0] : pSrc;
#         var["gSrc"] = pSrc;
#         var["gDst"] = pDst;
#         var["gLinearSampler"] = mpLinearSampler;
#         mpPostFXPass->execute(pRenderContext, uint3(resolution, 1));
#     }
# }

# void SimplePostFX::preparePostFX(RenderContext* pRenderContext, uint32_t width, uint32_t height)
# {
#     for (int res = 0; res < kNumLevels + 1; ++res)
#     {
#         ref<Texture>& pBuf = mpPyramid[res];
#         if (getBloomAmount() <= 0.f)
#         {
#             pBuf = nullptr;
#         }
#         else
#         {
#             uint32_t w = std::max(1u, width >> res);
#             uint32_t h = std::max(1u, height >> res);
#             if (!pBuf || pBuf->getWidth() != w || pBuf->getHeight() != h)
#             {
#                 pBuf = mpDevice->createTexture2D(
#                     w, h, ResourceFormat::RGBA16Float, 1, 1, nullptr, ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess
#                 );
#                 FALCOR_ASSERT(pBuf);
#             }
#         }
#     }
# }

# void SimplePostFX::renderUI(Gui::Widgets& widget)
# {
#     // Controls for output size.
#     // When output size requirements change, we'll trigger a graph recompile to update the render pass I/O sizes.
#     if (widget.dropdown("Output size", mOutputSizeSelection))
#         requestRecompile();
#     if (mOutputSizeSelection == RenderPassHelpers::IOSize::Fixed)
#     {
#         if (widget.var("Size in pixels", mFixedOutputSize, 32u, 16384u))
#             requestRecompile();
#     }

#     // PostFX options.
#     widget.checkbox("Enable post fx", mEnabled);
#     widget.slider("Wipe", mWipe, 0.f, 1.f);
#     if (auto group = widget.group("Lens FX", true))
#     {
#         group.slider("Bloom", mBloomAmount, 0.f, 1.f);
#         group.slider("Bloom Star", mStarAmount, 0.f, 1.f);
#         group.slider("Star Angle", mStarAngle, 0.f, 1.f, true);
#         group.slider("Vignette", mVignetteAmount, 0.f, 1.f);
#         group.slider("Chromatic Aberration", mChromaticAberrationAmount, 0.f, 1.f);
#         group.slider("Barrel Distortion", mBarrelDistortAmount, 0.f, 1.f);
#         if (group.button("reset this group"))
#         {
#             mBloomAmount = 0.f;
#             mStarAmount = 0.f;
#             mStarAngle = 0.1f;
#             mVignetteAmount = 0.f;
#             mChromaticAberrationAmount = 0.f;
#             mBarrelDistortAmount = 0.f;
#         }
#     }
#     if (auto group = widget.group("Saturation", true))
#     {
#         group.slider("Shadow Saturation", mSaturationCurve.x, 0.f, 2.f);
#         group.slider("Midtone Saturation", mSaturationCurve.y, 0.f, 2.f);
#         group.slider("Hilight Saturation", mSaturationCurve.z, 0.f, 2.f);
#         if (group.button("reset this group"))
#         {
#             mSaturationCurve = float3(1.f);
#         }
#     }
#     if (auto group = widget.group("Offset/Power/Scale (luma)", true))
#     {
#         group.slider("Luma Offset (Shadows)", mColorOffsetScalar, -1.f, 1.f);
#         group.slider("Luma Power (Midtones)", mColorPowerScalar, -1.f, 1.f);
#         group.slider("Luma Scale (Hilights)", mColorScaleScalar, -1.f, 1.f);
#         if (group.button("reset this group"))
#         {
#             mColorOffsetScalar = 0.f;
#             mColorPowerScalar = 0.f;
#             mColorScaleScalar = 0.f;
#         }
#     }
#     if (auto group = widget.group("Offset/Power/Scale (color)", true))
#     {
#         if (group.button("reset##1"))
#             mColorOffset = float3(0.5f, 0.5f, 0.5f);
#         group.rgbColor("Color Offset (Shadows)", mColorOffset, true);

#         if (group.button("reset##2"))
#             mColorPower = float3(0.5f, 0.5f, 0.5f);
#         group.rgbColor("Color Power (Midtones)", mColorPower, true);

#         if (group.button("reset##3"))
#             mColorScale = float3(0.5f, 0.5f, 0.5f);
#         group.rgbColor("Color Scale (Hilights)", mColorScale, true);

#         if (group.button("reset this group"))
#         {
#             mColorOffset = float3(0.5f, 0.5f, 0.5f);
#             mColorPower = float3(0.5f, 0.5f, 0.5f);
#             mColorScale = float3(0.5f, 0.5f, 0.5f);
#         }
#     }
# }