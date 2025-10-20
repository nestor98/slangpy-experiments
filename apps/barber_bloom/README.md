# Shadertoy + screenspace postprocessing
This is a test app that adds a bloom pass (and some other postFX from falcor that simplify something like the compositing nodes in blender) to the barber pole shader. This allows removing removing the per-pixel minimum distance array. That should run faster ?? 

Most heavy lifting is made by an adapted version of falcor's [simplePostFX](https://github.com/NVIDIAGameWorks/Falcor/blob/eb540f6748774680ce0039aaf3ac9279266ec521/Source/RenderPasses/SimplePostFX/SimplePostFX.cs.slang) slang file. 

Update: 
- the postFX works and is quite modular, check out the code to add it to other apps.
- It ofc looks a bit different (bc bloom is only screenspace while the old effect was 3d), but it also allows some other nice stuff like vignetting & color grading in a modular way
- TODO: measure performance vs original
