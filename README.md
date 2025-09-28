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

### Barber pole ([link](https://github.com/nestor98/slangpy-experiments/tree/master/apps/barber1))
- Ported my [shadertoy visualization](https://www.shadertoy.com/view/3fsGzr) of the barber pole effect. It's a good example for porting other shaders. There were some gotchas in the process, but with slang and python some things (such as inputs and managing state) are much less painful than in shadertoy. 
- It was made roughly following [this example](https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang) about wrapping Shadertoy from the slang repo, but removing the interface stuff... Didn't seem very useful and it was quite complicated
- Added a couple of sliders
- It would be interesting to think how we could make use of compute shaders here to speed things up
