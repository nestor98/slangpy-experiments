## Some random slang experiments

### Barber pole ()
- Ported my shadertoy visualization of the barber pole effect. It's a good example for porting other shaders. There were some gotchas in the process, but with slang and python some things (such as inputs and managing state) are much less painful than in shadertoy). 
- It was made roughly following [this example](https://github.com/shader-slang/slang/blob/master/examples/shader-toy/shader-toy.slang) from the slang repo, but removing the interface stuff... Didn't seem very useful
- Added a couple of sliders
- It would be interesting to think how we could make use of compute shaders here to speed things up
