## Some random slang experiments

### Installation
Make a venv and activate it (or conda). Install:
```
pip install -r requirements.txt
```
Run something: 
```
python apps/barber1/barberpole.py apps/barber1/shaders/barber.glsl
```

the shadertoy based apps preprocess glsl to slang. This may not work in all cases without tweaking the glsl or the generated slang, but is straightforward for simple shadertoys. 

You can also directly input the transpiled slang to avoid redoing the preprocessing:

```
python apps/barber1/barberpole.py apps/barber1/wrapped_shader.slang
```

### Barber pole ([src](https://github.com/nestor98/slangpy-experiments/tree/master/apps/barber1))
```
python apps/barber1/barberpole.py apps/barber1/wrapped_shader.slang
```
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/8dba403e-1dc0-4d67-81ef-794ecbf2a1f6" />

- Ported my [shadertoy visualization](https://www.shadertoy.com/view/3fsGzr) of the barber pole effect. It's a good example for porting other shaders. There were some gotchas in the process, but with slang and python some things (such as inputs and managing state) are much less painful than in shadertoy. 
- It was made roughly following [this example](https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang) about wrapping Shadertoy from the slang repo, but removing the interface stuff... Didn't seem very useful and it was quite complicated
- Added a couple of sliders
- It would be interesting to think how we could make use of compute shaders here to speed things up

### Barber pole with PostFX ([src](https://github.com/nestor98/slangpy-experiments/tree/master/apps/barber_bloom))
```
python apps/barber_bloom/barberpole.py 
```
<img width="1715" height="857" alt="Captura de pantalla 2025-09-29 233710" src="https://github.com/user-attachments/assets/4c47b7a8-e8db-4222-bcd6-8cb8a3e590c1" />
<img width="1920" height="1020" alt="Captura de pantalla 2025-09-29 235528" src="https://github.com/user-attachments/assets/b04fc53b-4b5c-4496-adb5-5957d32f1317" />


- Adapted that to run with a modified version of a Falcor postprocessing .slang file (see [the Falcor source](https://github.com/NVIDIAGameWorks/Falcor/blob/eb540f6748774680ce0039aaf3ac9279266ec521/Source/RenderPasses/SimplePostFX/SimplePostFX.cs.slang), by Simon Kallweit)
- Fixed some artifacts in the shader (NaNs and SDF false positives)
- Bloom can still be adapted to look more like the original, i guess

