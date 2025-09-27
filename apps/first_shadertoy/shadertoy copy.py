import sys
import time
from pathlib import Path
import numpy as np
import slangpy as spy

# -----------------------
# User shader wrapper
# -----------------------

USER_TEMPLATE = """
import Shadertoy;

struct UserShader : IShaderToyImageShader
{
    {user_code}

    static This getDefault() {{ return This(); }}
};
"""

def wrap_user_shader(src: str) -> str:
    return USER_TEMPLATE.format(user_code=src)


# -----------------------
# App
# -----------------------

import sys
import time
from pathlib import Path
import numpy as np
import slangpy as spy

EXAMPLE_DIR = Path(__file__).parent
USE_RAYTRACING_PIPELINE = False  # adjust if needed

# -----------------------
# ShaderToy wrapper
# -----------------------

SHADER_TEMPLATE = """
import Shadertoy;


struct UserShader : IShaderToyImageShader {{
    {user_code}

    static This getDefault() {{ return This(); }}
}};
"""

def wrap_user_shader(src: str) -> str:
    return SHADER_TEMPLATE.format(user_code=src)

# -----------------------
# App
# -----------------------

class App:
    def __init__(self, user_shader_path: str):
        super().__init__()
        self.window = spy.Window(width=1920, height=1080, title="ShaderToy", resizable=True)
        self.device = spy.Device(
            enable_debug_layers=False,
            compiler_options={
                "include_paths": [EXAMPLE_DIR],
                "defines": {"USE_RAYTRACING_PIPELINE": "1" if USE_RAYTRACING_PIPELINE else "0"},
            },
        )
        self.surface = self.device.create_surface(self.window)
        self.surface.configure(width=self.window.width, height=self.window.height, vsync=False)

        self.output_texture: spy.Texture = None
        self.render_texture: spy.Texture = None
        self.accum_texture: spy.Texture = None

        self.window.on_keyboard_event = self.on_keyboard_event
        self.window.on_mouse_event = self.on_mouse_event
        self.window.on_resize = self.on_resize

        # Load user shader
        src = Path(user_shader_path).read_text()
        full_src = wrap_user_shader(src)
        self.module = spy.Module.load_from_source(self.device, "UserShader", full_src)

        # Dummy 1x1 textures for channels
        self.channels = [
            self.device.create_texture(1, 1, np.zeros((1, 1, 4), dtype=np.uint8)) for _ in range(4)
        ]

        self.start_time = time.time()

    def on_keyboard_event(self, event: spy.KeyboardEvent):
        if event.type == spy.KeyboardEventType.key_press:
            if event.key == spy.KeyCode.escape:
                self.window.close()

    def on_mouse_event(self, event: spy.MouseEvent):
        pass  # extend as needed

    def on_resize(self, width: int, height: int):
        self.device.wait()
        if width > 0 and height > 0:
            self.surface.configure(width=width, height=height, vsync=False)
        else:
            self.surface.unconfigure()

    def main_loop(self):
        frame = 0
        timer = spy.Timer()

        while not self.window.should_close():
            dt = timer.elapsed_s()
            timer.reset()

            self.window.process_events()
            if not self.surface.config:
                continue

            surface_texture = self.surface.acquire_next_image()
            if not surface_texture:
                continue

            # Allocate output textures if needed
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
                self.render_texture = self.device.create_texture(
                    format=spy.Format.rgba32_float,
                    width=surface_texture.width,
                    height=surface_texture.height,
                    usage=spy.TextureUsage.shader_resource | spy.TextureUsage.unordered_access,
                )
                self.accum_texture = self.device.create_texture(
                    format=spy.Format.rgba32_float,
                    width=surface_texture.width,
                    height=surface_texture.height,
                    usage=spy.TextureUsage.shader_resource | spy.TextureUsage.unordered_access,
                )

            command_encoder = self.device.create_command_encoder()

            # Dispatch ShaderToy shader
            kwargs = {
                "iResolution": (float(self.window.width), float(self.window.height)),
                "iTime": float(time.time() - self.start_time),
                "iMouse": self.window.mouse,
                "iChannel0": self.channels[0],
                "iChannel1": self.channels[1],
                "iChannel2": self.channels[2],
                "iChannel3": self.channels[3],
            }

            self.module.dispatch(command_encoder, self.render_texture, **kwargs)

            # Copy to output texture
            command_encoder.blit(self.output_texture, self.render_texture)
            self.device.submit_command_buffer(command_encoder.finish())
            del surface_texture

            self.surface.present()
            frame += 1

# -----------------------
# CLI
# -----------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python slangpy_shadertoy_app.py path/to/shader.glsl")
        sys.exit(1)

    app = App(sys.argv[1])
    app.main_loop()
