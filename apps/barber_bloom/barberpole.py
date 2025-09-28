import re
import sys
import time
from pathlib import Path
import numpy as np
import slangpy as spy

BASE_DIR = Path(__file__).parent



def wrap_shadertoy_glsl(src: str, slang_file_1="shadertoy1.slang", slang_file_2="shadertoy2.slang") -> str:
    """Wraps GLSL code into Slang compute shader with ShaderToy uniforms."""
    # I cannot get the import to work, so doing it manually here
    shadertoy_defs_src = Path(__file__).parent / slang_file_1
    shadertoy_defs_src = shadertoy_defs_src.read_text()

    shadertoy_wrapped_shader =  Path(__file__).parent / slang_file_2
    shadertoy_wrapped_shader = shadertoy_wrapped_shader.read_text()
    shadertoy_wrapped_shader = shadertoy_wrapped_shader.replace("{USER_CODE}", src)

    source = shadertoy_defs_src + "\n" + shadertoy_wrapped_shader
    return source

# -----------------------
# Runtime App
# -----------------------
class ShaderToyUniforms:
    def __init__(self, resolution, time, mouse):
        self.iResolution = resolution  # float2
        self.iTime = time             # float
        self.iMouse = mouse           # float4


def glsl_src_to_slang(glsl_src, slang_file_1="shadertoy1.slang", 
                      slang_file_2="shadertoy2.slang",
                      output_slang="wrapped_shader.slang") -> str:
    """
    Converts GLSL source to Slang source by wrapping it with slang shadertoy defines.
    Returns the slang source, and also writes it to output_slang file.
    """
    slang_src = wrap_shadertoy_glsl(glsl_src)
    temp_slang_file = BASE_DIR / "wrapped_shader.slang"
    with temp_slang_file.open("w") as f:
        f.write(slang_src)
    return str(slang_src)

def glsl_file_to_slang(glsl_file, slang_file_1="shadertoy1.slang", 
                      slang_file_2="shadertoy2.slang",
                      output_slang="wrapped_shader.slang") -> str:
    """ Same as above, but reads the glsl from a file """
    glsl_src = Path(glsl_file).read_text()
    return glsl_src_to_slang(glsl_src, slang_file_1, slang_file_2, output_slang)


class App:
    def __init__(self, src_file=None, flip_mouse_y=True):
        super().__init__()
        self.window = spy.Window(width=1920, height=1080, title="ShaderToy", resizable=True)
        # self.device = spy.Device(enable_debug_layers=False, compiler_options={"include_paths": [BASE_DIR]})
                
        self.device = spy.create_device(
            include_paths=[
                Path(__file__).parent.absolute(),
            ],
            enable_debug_layers=False # Debug is a lot slower (more than 2x)
        )

        self.surface = self.device.create_surface(self.window)
        self.surface.configure(width=self.window.width, height=self.window.height)

        self.ui = spy.ui.Context(self.device)

        self.output_texture: spy.Texture = None

        self.logger = spy.Logger()

        # check if the source file is glsl or slang
        src_path = Path(src_file)
        if not src_path.exists():
            self.logger.error(f"File {src_file} does not exist")
            sys.exit(1)
        if src_path.suffix not in [".glsl", ".slang"]:
            self.logger.error(f"File {src_file} is not a .glsl or .slang file")
            sys.exit(1)
        if src_path.suffix == ".slang":
            self.logger.info("Using .slang file directly (must have the correct ShaderToyUniforms and entry point).")
            temp_slang_file = src_file
        else:
            # Load and wrap GLSL
            glsl_src = Path(src_file).read_text()
            slang_src = wrap_shadertoy_glsl(glsl_src)
            temp_slang_file = BASE_DIR / "wrapped_shader.slang"
            with temp_slang_file.open("w") as f:
                f.write(slang_src)
            self.logger.info(f"Wrapped .glsl file into {temp_slang_file}")

        self.program = self.device.load_program(str(temp_slang_file), ["compute_main"])
        
        self.kernel = self.device.create_compute_kernel(self.program)

        self.start_time = time.time()
        self.frame = 0

        # Dummy channel textures
        self.channels = [
            self.device.create_texture(
                format=spy.Format.rgba32_float,
                width=1, height=1, data=np.zeros((1, 1, 4), dtype=np.float32))
            for _ in range(4)
        ]

        self.R = spy.float2(self.window.width, self.window.height)
        self.mouse_pos = spy.float2(.5) * self.R 
        self.mouse_pos_shift = spy.float2(.5) * self.R 
        self.mouse_pos_ctrl = spy.float2(.5) * self.R 
        self.mouse_pos_alt = spy.float2(.5) * self.R 
        self.rmouse_pos = spy.float2(.5) * self.R 
        self.flip_mouse_y = flip_mouse_y

        self.mouse_down = False
        self.rmouse_down = False

        self.window.on_keyboard_event = self.on_keyboard_event
        self.window.on_mouse_event = self.on_mouse_event
        self.window.on_resize = self.on_resize

        self.playing=True
        self.fps_avg = 0.0
        
        self.setup_ui()

        # post fx:  
        # TODO: add the postfx module using the falcor slang file:
        # self.postfx_module = spy.Module("postfx", self.device)
        self.enable_postfx = False

    def on_keyboard_event(self, event: spy.KeyboardEvent):
        if self.ui.handle_keyboard_event(event):
            return
        if event.type == spy.KeyboardEventType.key_press:
            if event.key == spy.KeyCode.escape:
                self.window.close()
            if event.key == spy.KeyCode.f11:
                spy.Warning("Fullscreen toggle not implemented")
                print("I wish i could toggle fullscreen, dk how :(")
            #     self.window.toggle_fullscreen()
            if event.key == spy.KeyCode.f:
                print(f"Frame {self.frame} Time {time.time() - self.start_time:.2f}s FPS {self.frame / (time.time() - self.start_time):.2f}")
            if event.key == spy.KeyCode.f2:
                spy.tev.show_async(self.output_texture, "shadertoy_output")
            if event.key == spy.KeyCode.space:
                self.ui_window.set_visible(not self.ui_window.visible())
                

    def setup_ui(self):
        screen = self.ui.screen
        def toggle_settings():
            self.ui_window.show()
            ui_button.close()
        # ui_button = spy.ui.Button(label="?", callback=toggle_settings)
        # ui_button.tooltip = "Toggle Settings"

        self.ui_window = spy.ui.Window(screen, "Settings", size=spy.float2(150, 60))
        # window.close()
        window = self.ui_window

        self.fps_text = spy.ui.Text(window, "FPS: 0")

        def play():
            self.playing = not self.playing

        spy.ui.Button(window, "Play/Pause", callback=play)

        self.freq_mult = spy.ui.SliderFloat(window, "Frequency", value=1, min=0.001, max=5)
        self.twist_mult = spy.ui.SliderFloat(window, "Twist", value=1, min=0, max=3)
        # self.mouse_radius = spy.ui.SliderFloat(window, "Radius", value=100, min=0, max=1000)

        spy.ui.Button(window, "Toggle PostFX", callback=lambda: setattr(self, 'enable_postfx', not self.enable_postfx))

    def on_mouse_event(self, event: spy.MouseEvent):
        if self.ui.handle_mouse_event(event): # if mouse is over ui, it has priority
            return
        
        if self.flip_mouse_y:
            event.pos.y = self.window.height - event.pos.y

        if event.type == spy.MouseEventType.move:
            if event.has_modifier(spy.KeyModifier.shift):
                self.mouse_pos_shift = event.pos
            elif event.has_modifier(spy.KeyModifier.ctrl):
                self.mouse_pos_ctrl = event.pos
            elif event.has_modifier(spy.KeyModifier.alt):
                self.mouse_pos_alt = event.pos
            elif self.mouse_down:
                self.mouse_pos = event.pos
            elif self.rmouse_down:
                self.rmouse_pos = event.pos
        elif event.type == spy.MouseEventType.button_down:
            if event.button == spy.MouseButton.left:
                self.mouse_down = True
            elif event.button == spy.MouseButton.right:
                self.rmouse_down = True
        elif event.type == spy.MouseEventType.button_up:
            if event.button == spy.MouseButton.left:
                self.mouse_down = False
            elif event.button == spy.MouseButton.right:
                self.rmouse_down = False

    def on_resize(self, width: int, height: int):
        self.device.wait()
        if width > 0 and height > 0:
            self.surface.configure(width=width, height=height)
        else:
            self.surface.unconfigure()

    def run(self):
        timer = spy.Timer()
        sim_time = 0.0

        while not self.window.should_close():
            self.window.process_events()
            self.ui.process_events()

            elapsed = timer.elapsed_s()
            timer.reset()
            if self.playing:
                sim_time += elapsed
                
            self.fps_avg = 0.95 * self.fps_avg + 0.05 * (1.0 / elapsed)
            self.fps_text.text = f"FPS: {self.fps_avg:.2f}"

            if not self.surface.config:
                continue

            surface_texture = self.surface.acquire_next_image()
            if not surface_texture:
                continue

            if (
                self.output_texture is None
                or self.output_texture.width != surface_texture.width
                or self.output_texture.height != surface_texture.height
            ):
                self.output_texture = self.device.create_texture(
                    format=spy.Format.rgba32_float,
                    width=surface_texture.width,
                    height=surface_texture.height,
                    usage=spy.TextureUsage.shader_resource | spy.TextureUsage.unordered_access,

                )

            command_encoder = self.device.create_command_encoder()

            self.kernel.dispatch(
                thread_count=[self.output_texture.width, self.output_texture.height, 1],
                vars={
                    "g_output": self.output_texture,
                    "ShaderToyUniforms" : { # constant buffer passed as dict
                        "iResolution":(float(self.window.width), float(self.window.height)),
                        "iTime": float(sim_time)+10,
                        "iMouse": (self.mouse_pos.x, self.mouse_pos.y, 1.0 if self.mouse_down else 0.0, 0.0),
                        "iMouseRight": (self.rmouse_pos.x, self.rmouse_pos.y, 1.0 if self.rmouse_down else 0.0, 0.0),
                        "iMouseShift": (self.mouse_pos_shift.x, self.mouse_pos_shift.y, 0.0, 0.0),
                        "iMouseCtrl": (self.mouse_pos_ctrl.x, self.mouse_pos_ctrl.y, 0.0, 0.0),
                        "iMouseAlt": (self.mouse_pos_alt.x, self.mouse_pos_alt.y, 0.0, 0.0),
                        # specific for this shader:
                        "g_freq_mult": self.freq_mult.value,
                        "g_twist_mult": self.twist_mult.value,
                    },
                    # "iChannel0": self.channels[0],
                    # "iChannel1": self.channels[1],
                    # "iChannel2": self.channels[2],
                    # "iChannel3": self.channels[3],
                },
                command_encoder=command_encoder,
            )

            command_encoder.blit(surface_texture, self.output_texture)
            # ui:
            self.ui.new_frame(surface_texture.width, surface_texture.height)
            self.ui.render(surface_texture, command_encoder)
            # ................
            self.device.submit_command_buffer(command_encoder.finish())
            del surface_texture

            self.surface.present()
            self.frame += 1


if __name__ == "__main__":
    slang_file = Path(__file__).parent / "barber.slang"
    if len(sys.argv) > 1:
        slang_file = Path(sys.argv[1])

    app = App(str(slang_file))
    app.run()
