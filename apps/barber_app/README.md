
The version in barber1 is lower level, using devices and so on directly. It was made following examples such as https://github.com/shader-slang/slangpy/tree/main/examples/pathtracer

I later realized that the slangpy-samples examples seem to use a "newer" API, because slangpy slowly became higher level. For example, in experiments such as https://github.com/shader-slang/slangpy-samples/tree/main/experiments/brdf, they wrap code in an App and abstract a bit more devices and use more sp.Module 

I will try to do something like this here

