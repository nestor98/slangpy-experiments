import re
import sys
import time
from pathlib import Path
import numpy as np
import slangpy as spy

BASE_DIR = Path(__file__).parent



def wrap_shadertoy_glsl(src: str) -> str:
    """Wraps GLSL code into Slang compute shader with ShaderToy uniforms."""
    # I cannot get the import to work, so doing it manually here
    shadertoy_defs_src = Path(__file__).parent / "shadertoy1.slang"
    shadertoy_defs_src = shadertoy_defs_src.read_text()

    shadertoy_wrapped_shader =  Path(__file__).parent / "shadertoy2.slang"
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


class App:
    def __init__(self, glsl_path: str):
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

        self.output_texture: spy.Texture = None

        # Load and wrap GLSL
        glsl_src = Path(glsl_path).read_text()
        slang_src = wrap_shadertoy_glsl(glsl_src)
        temp_slang_file = BASE_DIR / "wrapped_shader.slang"
        with temp_slang_file.open("w") as f:
            f.write(slang_src)


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

        self.mouse_pos = spy.float2()
        self.mouse_down = False

        self.window.on_keyboard_event = self.on_keyboard_event
        self.window.on_mouse_event = self.on_mouse_event
        self.window.on_resize = self.on_resize

    def on_keyboard_event(self, event: spy.KeyboardEvent):
        if event.type == spy.KeyboardEventType.key_press:
            if event.key == spy.KeyCode.escape:
                self.window.close()
            if event.key == spy.KeyCode.f11:
                self.window.toggle_fullscreen()
            if event.key == spy.KeyCode.f:
                print(f"Frame {self.frame} Time {time.time() - self.start_time:.2f}s FPS {self.frame / (time.time() - self.start_time):.2f}")
            if event.key == spy.KeyCode.f2:
                spy.tev.show_async(self.output_texture, "shadertoy_output")

    def on_mouse_event(self, event: spy.MouseEvent):
        if event.type == spy.MouseEventType.move:
            self.mouse_pos = event.pos
        elif event.type == spy.MouseEventType.button_down:
            if event.button == spy.MouseButton.left:
                self.mouse_down = True
        elif event.type == spy.MouseEventType.button_up:
            if event.button == spy.MouseButton.left:
                self.mouse_down = False

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
            elapsed = timer.elapsed_s()
            timer.reset()
            sim_time += elapsed

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
                        "iTime": float(sim_time),
                        "iMouse": (self.mouse_pos.x, self.mouse_pos.y, 1.0 if self.mouse_down else 0.0, 0.0),
                    },
                    # "iChannel0": self.channels[0],
                    # "iChannel1": self.channels[1],
                    # "iChannel2": self.channels[2],
                    # "iChannel3": self.channels[3],
                },
                command_encoder=command_encoder,
            )

            command_encoder.blit(surface_texture, self.output_texture)
            self.device.submit_command_buffer(command_encoder.finish())
            del surface_texture

            self.surface.present()
            self.frame += 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python slangpy_shadertoy_app.py path/to/shadertoy.glsl")
        sys.exit(1)

    app = App(sys.argv[1])
    app.run()
