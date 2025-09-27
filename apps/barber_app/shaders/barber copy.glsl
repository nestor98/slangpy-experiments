// source: https://www.shadertoy.com/view/flXBzB
#define PI 3.141592653

// old defines:
#define MOUSE_MODE 3 
#define VIEW 3
#define REPEAT_DOMAIN 0


// Some defines had to be moved to shadertoy2.slang for language limitations >:(


//vec2 iCam; 




// original: https://www.shadertoy.com/view/7djSzK
// See: https://www.desmos.com/calculator/hvqys18zux
// f: frequency
// a: amplitude
float sdSine(in vec2 p, in float f, in float a) {
    f *= 3.14159265359 * a; p /= a; // Modify to handle varying amplitude
    float r = 3.14159265359 / f, h = 0.5 * r, ff = f * f;
    p = vec2(mod(p.x + h, r) - h, p.y * sign(r - mod(p.x + h, 2.0 * r))); // Remap

    // Get closest on linear approximation
    float t = clamp((0.818309886184 * f * p.y + p.x) / (0.669631069826 * ff + 1.0), -h, h);

    // Three iterations of Newton-Raphson
    for (int n=0; n < 3; n++) {
        float k = t * f, c = cos(k), s = sin(k);
        t -= ((s - p.y) * c * f + t - p.x) / ((c * c - s * s + s * p.y) * ff + 1.0);
    }

    return length(p - vec2(t, sin(t * f))) * a;
}

// -------------------- CIE XYZ 2006 ----------------------------

vec3 xyz2srgb(vec3 xyz) {
    // From https://github.com/tobspr/GLSL-Color-Spaces
    const mat3 XYZ_2_RGB = (mat3(
         3.2404542,-0.9692660, 0.0556434,
        -1.5371385, 1.8760108,-0.2040259,
        -0.4985314, 0.0415560, 1.0572252
    ));
    return XYZ_2_RGB * xyz;
}


#ifndef OK_LAB
// t should be in 0-1. Returns the srgb color corresponding to that wl
// according to cie 2006
// super overkill for this app but why not
vec3 rgb_wl(float t) {
    const int ignore = 0; // too much red
    //t /= float(CIEXYZ_SAMPLES+ignore)/float(CIEXYZ_SAMPLES);
    int idx = int(t * float(CIEXYZ_SAMPLES-1));
    float interp = fract(t * float(CIEXYZ_SAMPLES-1));
    vec3 rgb = xyz2srgb(mix(ciexyz[idx], ciexyz[idx+1], interp));
    return rgb;
}


#else
//

// oklab = (lightness, red_greenness, blue_yelowness)
vec3 oklab2lrgb(vec3 oklab) {
    vec3 lms = oklab * mat3(1,  0.3963377774,  0.2158037573,
                            1, -0.1055613458, -0.0638541728,
                            1, -0.0894841775, -1.2914855480);
    lms *= lms * lms;
    return lms * mat3( 4.0767416621, -3.3077115913,  0.2309699292, 
                      -1.2684380046,  2.6097574011, -0.3413193965, 
                      -0.0041960863, -0.7034186147,  1.7076147010);
}

vec3 lrgb2oklab(vec3 lrgb) {
    vec3 lms = lrgb * mat3(0.4121656120, 0.5362752080, 0.0514575653,
                           0.2118591070, 0.6807189584, 0.1074065790,
                           0.0883097947, 0.2818474174, 0.6302613616);
    return pow(lms, vec3(1.0 / 3.0)) * mat3(0.2104542553,  0.7936177850, -0.0040720468,
                                            1.9779984951, -2.4285922050,  0.4505937099,
                                            0.0259040371,  0.7827717662, -0.8086757660);
}


// t should be in 0-1. Returns the srgb color corresponding to that wl
// according to cie 2006
// super overkill for this app but why not
vec3 rgb_wl(float t) {
    const int ignore = 0; // too much red
    //t /= float(CIEXYZ_SAMPLES+ignore)/float(CIEXYZ_SAMPLES);
    int idx = int(t * float(CIEXYZ_SAMPLES-1));
    float interp = fract(t * float(CIEXYZ_SAMPLES-1));
    vec3 rgb = xyz2srgb(mix(ciexyz[idx], ciexyz[idx+1], interp));
    return rgb;
}

#endif

/// Intersect plane (https://iquilezles.org/articles/intersectors/)
// plane degined by p (p.xyz must be normalized)
float plaIntersect( in vec3 ro, in vec3 rd, in vec4 p )
{
    return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}

// Pristine grid: https://www.shadertoy.com/view/mdVfWw
// version with explicit gradients for use with raycast shaders like this one
float pristineGrid( in vec2 uv, in vec2 ddx, in vec2 ddy, vec2 lineWidth)
{
    vec2 uvDeriv = vec2(length(vec2(ddx.x, ddy.x)), length(vec2(ddx.y, ddy.y)));
    bvec2 invertLine = bvec2(lineWidth.x > 0.5, lineWidth.y > 0.5);
    vec2 targetWidth = vec2(
      invertLine.x ? 1.0 - lineWidth.x : lineWidth.x,
      invertLine.y ? 1.0 - lineWidth.y : lineWidth.y
      );
    vec2 drawWidth = clamp(targetWidth, uvDeriv, vec2(0.5));
    vec2 lineAA = uvDeriv * 1.5;
    vec2 gridUV = abs(fract(uv) * 2.0 - 1.0);
    gridUV.x = invertLine.x ? gridUV.x : 1.0 - gridUV.x;
    gridUV.y = invertLine.y ? gridUV.y : 1.0 - gridUV.y;
    vec2 grid2 = smoothstep(drawWidth + lineAA, drawWidth - lineAA, gridUV);

    grid2 *= clamp(targetWidth / drawWidth, 0.0, 1.0);
    grid2 = mix(grid2, targetWidth, clamp(uvDeriv * 2.0 - 1.0, 0.0, 1.0));
    grid2.x = invertLine.x ? 1.0 - grid2.x : grid2.x;
    grid2.y = invertLine.y ? 1.0 - grid2.y : grid2.y;
    return mix(grid2.x, 1.0, grid2.y);
}

// By iq: https://iquilezles.org/articles/intersectors/
// cylinder defined by extremes a and b, and radious ra
vec4 cylIntersect( in vec3 ro, in vec3 rd, in vec3 a, in vec3 b, float ra )
{
    vec3  ba = b  - a;
    vec3  oc = ro - a;
    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoc = dot(ba,oc);
    float k2 = baba            - bard*bard;
    float k1 = baba*dot(oc,rd) - baoc*bard;
    float k0 = baba*dot(oc,oc) - baoc*baoc - ra*ra*baba;
    float h = k1*k1 - k2*k0;
    if( h<0.0 ) return vec4(-1.0);//no intersection
    h = sqrt(h);
    float t = (-k1-h)/k2;
    // body
    float y = baoc + t*bard;
    if( y>0.0 && y<baba ) return vec4( t, (oc+t*rd - ba*y/baba)/ra );
    // caps
    t = ( ((y<0.0) ? 0.0 : baba) - baoc)/bard;
    if( abs(k1+k2*t)<h )
    {
        return vec4( t, ba*sign(y)/sqrt(baba) );
    }
    return vec4(-1.0);//no intersection
}

// by iq: https://iquilezles.org/articles/distfunctions/
float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

// by knarkowicz:
// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
vec3 aces(vec3 col) {
    col *= 0.6f;
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((col*(a*col+b))/(col*(c*col+d)+e), 0.0f, 1.0f);
}

/* 
// TO DO: i think using something like iq's soft shadow for the glow would remove
// some of the banding of my naive method 
float softshadow( in vec3 ro, in vec3 rd, float mint, float maxt, float w )
{
    float res = 1.0;
    float t = mint;
    for( int i=0; i<256 && t<maxt; i++ )
    {
        float h = map(ro + t*rd);
        res = min( res, h/(w*t) );
        t += clamp(h, 0.005, 0.50);
        if( res<-1.0 || t>maxt ) break;
    }
    res = max(res,-1.0);
    return 0.25*(1.0+res)*(1.0+res)*(2.0-res);
}*/

// Common above

// MAin below
// Barber Pole effect

// Controls:
// Number keys change the mode
// - Press 1: mouse rotates camera
// - Press 2: mouse rotates second polarization filter
// - Press 3: mouse pans camera along the tube  (default)
// - Press 4: mouse moves camera closer or farther


// Explanation:
// This is a visualization of the barber pole polarization effect.
// Light of different wavelengths goes through a water and sugar solution.
// Sugar twists the plane of polarization (the direction these waves oscillate)
// along the tube. The rate of this twist depends on the wavelength of each wave.
// Therefore, using a linear polarizer at each end of the tube, white light
// turns different colors depending on the angles of the polarizers 
// (and the length of the tube)
// This is because the intensity of each wavelength is reduced by the squared cosine 
// of its angle with respect to the second filter - see here how the waves that are more aligned 
// to the second filter end up with a higher amplitudde

// > See https://www.youtube.com/watch?v=QCX62YJCmGk for a much more detailed explanation

// ofc most constants here (twisting rate, frequency ratios...) are made up and not physically based at all
// by default only 3 wavelengths are shown for clarity, in reality the spectrum is continuous

// Acks:
// Started From template "RayMarching starting point" by Martijn Steinrucken
// The sine wave SDF is from oneshade (https://www.shadertoy.com/view/7djSzK)
// ^ I tried to make one then quickly realized it is not as simple as it seems
// The twisting of the plane and many other functions are from iq (https://iquilezles.org/articles/distfunctions/)
// the beautiful grids are from https://www.shadertoy.com/view/mdVfWw
// -------------------------------------------------


// add wls here, with wl_idx in 0-1 indicating from blue (380) to red (700 something) 
#define N_WLS 3
// ^ Try modifying this, but Don't go higher than the array in line 43

//#define NO_TWIST
// ^uncomment to remove sugar, aka stop the waves from twisting



// this array has to have at least N_WLS wls
float wl_idx[11] = float[11](
        0.01, // blue 
        0.99, // red
        0.33, // green
        .5,   // yellow
        0.6,   // orange
        0.1,
        0.2,
        0.4,
        0.7,
        0.65,
        0.85 // red
        );    
        
// -------------------------------------------------
// drawing options:
#define FREQ_MULT 1.5
#ifdef NO_TWIST
#define TWIST_MULT 0.
#else
#define TWIST_MULT 1.
#endif

// Draw dots in the beginning of the wave:
#define DOTS 
#define CYLINDER

// -------------------------------------------------
// Marching parameters:
// ideally, STEP_DECREASE should be 1. (faster marching)
// this can only be done if sdfs are not distorted 
// (my sdWave thing is a little bc of the twisting)
// https://iquilezles.org/articles/distfunctions/
#define STEP_DECREASE .91
#define MAX_STEPS 100
#define MAX_DIST 20.
#define SURF_DIST .01
#define TAU 6.283185
    

#define TIME iTime


float dot2(vec3 a) {return dot(a,a);} 
mat2 Rot(float a) {
    float s=sin(a), c=cos(a);
    return mat2(c, -s, s, c);
}

float sdBox(vec3 p, vec3 s) {
    p = abs(p)-s;
	return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.);
}
// twists around y
vec3 opTwistX( in vec3 p, float k )
{
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(p.x, m*p.yz);
    return q;
}

// Wave along X axis
// n and h define the plane of the wave (normal and height along that normal)
// freq and amp, the freq and amplitude of the wave
float sdWave(vec3 p, float h, vec3 n, float freq, float amp) {
    float dPlane = dot(p,n) + h; 
    // The distance to the plane in which the sine lives
    // if undistorted, then we can compute the distance of that closest point to the sine 
    // and use pythagoras to find the full 3D distance..
    vec3 pPlane = p - n * dPlane; // point in plane
    
    vec2 uv = pPlane.xy + vec2(TIME,0.); // animate
    
    float dSine = sdSine(uv, freq, amp); // 2D nice dist to sine
    
    dSine = sqrt(dPlane*dPlane+dSine*dSine) - .02; // pythagoras, -.02 adds a bit of thickness
    dSine = max(dSine, p.x-x_range.x); // cut off an end
    dSine = max(dSine, -p.x+x_range.y); // cut off an end
    
#ifdef DOTS
    float y = amp*sin(freq * PI * TIME); 
    dSine = min(dSine, // union of dSine (the wave) and
                length(p-vec3(x_range.x, y, 0.))-.08 // a dot
            );
    // another dot at the end of the wave segment (looks better without this i think):
    //float phi = (x_range.y-x_range.x) * freq*PI; // phase shift depends on distance
    //y = amp*sin(freq * PI * TIME - phi);
    //dSine = min(dSine, // union of dSine (the wave) and
    //            length(p-vec3(x_range.y, y, 0.))-.08 // a dot
    //        ); // uncomment to add back second set of dots
            
#endif
    return dSine; 
}
 

// x is distance
// y is ID (0. for blue, 1. for red, etc)
// [mutating] // bc it changes minDists
vec2 GetDist(vec3 p, in float filter_theta, inout float minDists[N_WLS]) {
    float d = 1000.;//
    float id =0.;
    
    for (int i=0;i<N_WLS;i++) {
        float t = wl_idx[i]; // 0 blue 1 red
        float twist = TWIST_MULT*mix(min_twist, max_twist, t);
        float freq  = FREQ_MULT*mix(max_f, min_f, t);
        vec3 q = p;
#if REPEAT_DOMAIN == 1
        float s = 4.;
        q.yz = q.yz - s*round(q.yz/s);
#endif
        q = opTwistX(q, twist);
        float d2 =sdWave(q, 0., vec3(0.,0.,1.), freq, 1.);
        // record the minimum dist to this wave along this ray, for the glow effect:
        minDists[i] = min(minDists[i], d2); 
        if (d2<d) {
            d= d2;
            id = t;
        }
#if VIEW >= 2
        // Compute the waves coming out of the second filter (linear pol)
        // same as twist function, but using only the pos of the second filter instead of the x pos:
        if (filter_theta < -PI/2.) filter_theta += 2.*PI;
        float s = sin(filter_theta-PI/2.), 
              c = cos(filter_theta-PI/2.);   
              
        mat2  m = mat2(c,-s,s,c);
        // amplitude is the cosine of the angle between the twist and the filter orientation
        float angle = twist*x_range.y - filter_theta+PI/2.;
        float amp = cos(angle) * cos(angle);
        amp = max(0.000001,amp);//Remove singularity when amp=0 (perfectly perpendicular filters)
        
        //float amp = abs(sin(-twist*x_range.y+filter_theta));
        q = vec3(p.x-x_range.y, m*p.yz);
        d2 = sdWave(q, 0., vec3(0.,0.,1.), freq, amp);
        minDists[i] = min(minDists[i], d2); 

        if (d2<d) {
            d= d2;
            id = t;
        }
#endif
    }
    return vec2(d,id);
}

// x is distance from ro, 
// y is minimum distance of all of the steps,
// z is id of closest
// [mutating]
vec3 intersect( in vec3 ro, in vec3 rd, in float filter_theta, inout float minDists[N_WLS] )
{
    float id = -1.;
	const float maxd = MAX_DIST;
	float h = 1.0;
    float t = 0.0;
    // opt: march up to the bounding box
    float mind = MAX_DIST;
    
    // right wall plane at x_range.y-4.
    vec3 p = ro;
    
    for( int i=0; i<MAX_STEPS; i++ )
    {
        if( h<SURF_DIST || t>maxd || p.x < x_range.y-6. ) break;
        p = ro+rd*t;
	    vec2 dist =  GetDist( p,filter_theta, minDists );
        h = dist.x; id = dist.y;
        t += h * STEP_DECREASE;
        mind = min(abs(h), mind);
    }

    if( t>maxd || p.x < x_range.y-6.) t=-1.0;
	
    return vec3(t, mind, id);
}

// warning: don't use this with the glowy global array

/*
[mutating]
vec3 GetNormal( in vec3 pos )
{
    vec3 eps = vec3(0.01,0.0,0.0);

	return normalize( vec3(
           GetDist(pos+eps.xyy).x - GetDist(pos-eps.xyy).x,
           GetDist(pos+eps.yxy).x - GetDist(pos-eps.yxy).x,
           GetDist(pos+eps.yyx).x - GetDist(pos-eps.yyx).x ) );
}*/


vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 
        f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = f*z,
        i = c + uv.x*r + uv.y*u;
    return normalize(i);
}


// Reads from one of our variables ([0:4], 0). Just for interaction
vec4 read(ivec2 id) {
    return .5;// TEMP texelFetch(iChannel0, id, 0);
}


//[mutating]
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    

    fragColor = vec4(iTime, iTime*2., iTime/3., 1.0);
    return;
        
    float minDists[N_WLS]; // for cheap glowy thing

    vec2 uv = (fragCoord-.5*iResolution.xy)/iResolution.y;
	vec2 m = read(ivec2(1,0)).xy;//iMouse.xy/iResolution.xy;
    // initialize glowy vector:
    for (int i=0;i<N_WLS;i++) { minDists[i] = MAX_DIST; }
    
    vec3 ro = vec3(0, 3, -3);

#if MOUSE_MODE == 1
    m = vec2(.5);
#elif MOUSE_MODE == 2
    m = vec2(.10,-.40);
#endif
    
    // ----------------------
    // Camera movement:
#if MOUSE_MODE == 3
    // move cam back or front
    ro *= 3. * (max(read(ivec2(4,0)).y, .01)); 
#endif

    ro.yz = ro.yz * Rot(-m.y*PI+1.);
    ro.xz = ro.xz * Rot(-m.x*TAU);

    vec3 rd = GetRayDir(uv, ro, vec3(0,0.,0), 1.);
    vec3 col = vec3(0);
    
    ro.x -= 2.; // move cam to right a bit
    
    #if VIEW >= 1
    ro.x -= -x_range.y-1.;//28.;
    #endif
    #if VIEW == 3
    ro.x -= 1.8; 
    ro.x -= mix(-11.8, 2., read(ivec2(3,0)).x);
    ro.x = max(ro.x, x_range.y-5.9); // keep the cam bounded before final wall
    #endif
    // ----------------------
    // Filter rotation:
        
    float filter_theta=3.*PI/4.;
    
    mat2 filter_R = mat2 (0.,-1.,1,0.);
#if MOUSE_MODE == 2
    filter_theta = atan(m.y-.5, m.x-.5);
#elif MOUSE_MODE == 3
    vec2 th = read(ivec2(2,0)).xy;
    filter_theta = -atan(th.y-.5, th.x-.5);
#endif
    float s = sin(filter_theta), 
          c = cos(filter_theta);
    filter_R = mat2(c,-s,s,c);
    // ----------------------
    
    vec3 intersection = intersect(ro, rd, filter_theta, minDists);
    float d = intersection.x; // dist origin->surface (-1 if no intersect)
    float closest_d = intersection.y; // for missing rays, how close did they get?
    float id = intersection.z; // id of closest thing
    
    // Shade the waves:
    vec3 p = ro+rd*MAX_DIST; // p of intersection
    float t = max(0., closest_d); // closest distance between ray & surface
    if (d>0.) { // did actually intersect
        t = 0.;
        p = ro+rd*d;
        
    }
    col = .9*vec3(.9,.9,.99) 
          * exp(-68.*t); // this is to shade even grazing rays (kind of AA)   
   
 #define DRAW_FILTERS
 #ifdef DRAW_FILTERS
    { // raytrace tube ends
        float tPlane = plaIntersect(ro,rd,vec4(1., 0., 0., -x_range.x));
        
        
        if (tPlane < MAX_DIST*5. && tPlane>0.) {
            vec3 pPlane=ro+rd*tPlane;
            vec2 ddx = dFdx(pPlane.yz), ddy = dFdy(pPlane.yz);
            
            float grid = pristineGrid( pPlane.yz, ddx, ddy, vec2(0.,0.01));
            
            // "lighting" for the grid:
            grid *= 30.*exp(-.56*length(pPlane.yz)); 
            
            col += .1*vec3(1.,1.,.8)*grid*(exp(-.02*tPlane)); 
        }
    }
    { // Second grid
        float tPlane = plaIntersect(ro,rd,vec4(1., 0., 0., -x_range.y));
        
        if (tPlane < MAX_DIST*5. && tPlane>0.){
            vec3 pPlane=ro+rd*tPlane;
            pPlane.yz = filter_R * pPlane.yz; // rotate the second filter
            
            vec2 ddx = dFdx(pPlane.yz), ddy = dFdy(pPlane.yz);
            
            float grid = pristineGrid( pPlane.yz, ddx, ddy, vec2(0.01,0.));
            
            // "lighting" for the grid:
            grid *= 30.*exp(-.56*length(pPlane.yz)); 
            
            col += .1*vec3(1.,1.,.8)*grid*(exp(-.02*tPlane)); 
        }
        
    }
 #endif
 #ifdef CYLINDER
 {
     vec4 cyl = cylIntersect( ro, rd, vec3(x_range.x,0.,0.), 
                         vec3(x_range.y,0.,0.), 1.6);
     
     float t = cyl.x;
     if (t>0.) {
         vec3 l = vec3(.7,1.,1.); // bluish
         vec3 n = normalize(cyl.yzw);
         #if 1
         l *= pow(1.-dot(n.yz, normalize(-rd.yz)),5.);//highlight the edge
         l *= 1.*smoothstep(.0,.3,dot(n.yz,normalize(-rd.yz))); // AA the edge 
         #else
         // also accounts for x normal (i think the other looks better)
         l *= pow(1.-dot(n, normalize(-rd)),11.);//highlight the edge
         l *= 3.*smoothstep(.0,.3,dot(n,normalize(-rd))); // AA the edge 
         #endif
         
         col += l;
     }
 }
 #endif
 
    // cheap glow using an imaginary extinction and the minimum distance found along each ray
    for (int i=0; i<N_WLS; i++){ // each wl has a color, and a min distance to this ray
        vec3 glow = normalize(rgb_wl(wl_idx[i]));
        float t = max(minDists[i],.0);
        col += .6*glow * exp(-15. * t);
        
    } 
 #if VIEW >= 2
 
    { // raytrace final wall
        float tPlane = plaIntersect(ro,rd,vec4(-1., 0., 0., x_range.y-6.));
        if (tPlane < MAX_DIST*5. && tPlane>0.) {
            vec3 pPlane=ro+rd*tPlane;
            
            // Get final color (avg of wls.)
            vec3 c = vec3(0.);
            
            vec2 ddx = dFdx(pPlane.yz), ddy = dFdy(pPlane.yz);
            float grid = pristineGrid( pPlane.yz, ddx, ddy, vec2(0.01,0.01));
            c += grid;

            if (filter_theta < -PI/2.) filter_theta += 2.*PI;
            for (int i=0;i<N_WLS;i++) {

                float t = wl_idx[i]; // 0 blue 1 red
                float twist = TWIST_MULT*mix(min_twist, max_twist, t);
                // amplitude is the cosine of the angle between the twist and the filter orientation
                float angle = twist*x_range.y - filter_theta+PI/2.;
                float amp = cos(angle) * cos(angle);
                
                c += amp*normalize(rgb_wl(wl_idx[i]));
            }
            c /= float(N_WLS);
            c = clamp(c,0.,1.);
            col += 3.*c * exp(-.70*length(pPlane.zy));
        }
            
    }
 #endif
 
    //col =1.1* col/(1.+col); // reinhard (dont like it here) 
    col = aces(1.*col);
    col = pow(col, vec3(1./2.2)); // gamma
    fragColor = vec4(col,1.0);
}