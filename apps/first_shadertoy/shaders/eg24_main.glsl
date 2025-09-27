
// -----------------------------------------------------------
// COMMON https://www.shadertoy.com/view/4cVGD3

// Data (spectra [Williamson&Hollins22: https://doi.org/10.1364/AO.470464] & sensitivity) 
// and utility functions, many many of which are the amazing work of iq

// The Jerlov water types we use here are useful for modeling a variety of the
// natural waters in the world
// For example, it seems like Weta also used them as a base for their waters in Avatar 2 
// [https://dl.acm.org/doi/abs/10.1145/3587423.3595519]
// -----------------------------------------------------------------------------

#define PI 3.141592

// -----------------------------------------------------------------------------
// Water optical properties, from Measured IOPs of Jerlov water types [Williamson and Hollins 22].
// Here I gathered 20 (linearly interpolated) wavelength samples, but the rest can be found here:
// https://doi.org/10.6084/m9.figshare.20290782

// Here, there are some examples of water optical properties from the jerlov water types
// Each vec3 contains: (scattering b, extinction c, Kd)
// (Where Kd is the diffuse downwelling attenuation coefficient)
// Wavelengths: linspace(400,700,20)

// Number of wls in our spectra arrays (dont modify, you can modify <samples> in Image tab):
#define N_WLS 20
// ----------------- Oceanic:
// Jerlov IB
const vec3 water_b_c_kd_IB[N_WLS] = vec3[](
    vec3(0.144,0.189,0.051),vec3(0.141,0.184,0.045),vec3(0.138,0.180,0.040),vec3(0.135,0.176,0.037),vec3(0.133,0.171,0.034),vec3(0.130,0.167,0.034),vec3(0.129,0.167,0.040),vec3(0.127,0.173,0.047),vec3(0.124,0.176,0.055),vec3(0.123,0.182,0.066),vec3(0.121,0.193,0.080),vec3(0.120,0.208,0.098),vec3(0.119,0.257,0.184),vec3(0.117,0.369,0.260),vec3(0.116,0.402,0.304),vec3(0.115,0.425,0.343),vec3(0.114,0.461,0.381),vec3(0.113,0.525,0.419),vec3(0.112,0.565,0.488),vec3(0.111,0.684,0.580)
);
// Jerlov II
const vec3 water_b_c_kd_II[N_WLS] = vec3[](
    vec3(0.209,0.283,0.096),vec3(0.204,0.275,0.087),vec3(0.200,0.267,0.078),vec3(0.196,0.259,0.069),vec3(0.192,0.250,0.065),vec3(0.189,0.243,0.063),vec3(0.186,0.239,0.068),vec3(0.183,0.240,0.073),vec3(0.180,0.241,0.077),vec3(0.178,0.244,0.085),vec3(0.175,0.253,0.097),vec3(0.173,0.267,0.114),vec3(0.171,0.314,0.199),vec3(0.169,0.424,0.276),vec3(0.167,0.456,0.323),vec3(0.165,0.479,0.366),vec3(0.164,0.514,0.407),vec3(0.162,0.578,0.448),vec3(0.161,0.618,0.518),vec3(0.159,0.734,0.610)
);
// Jerlov III
const vec3 water_b_c_kd_III[N_WLS] = vec3[](
    vec3(0.324,0.454,0.185),vec3(0.319,0.447,0.169),vec3(0.315,0.439,0.153),vec3(0.310,0.424,0.138),vec3(0.306,0.408,0.125),vec3(0.302,0.394,0.116),vec3(0.298,0.383,0.115),vec3(0.295,0.378,0.115),vec3(0.292,0.374,0.116),vec3(0.289,0.372,0.119),vec3(0.285,0.377,0.129),vec3(0.282,0.387,0.147),vec3(0.280,0.434,0.233),vec3(0.277,0.543,0.312),vec3(0.275,0.576,0.362),vec3(0.273,0.598,0.408),vec3(0.270,0.635,0.453),vec3(0.268,0.706,0.500),vec3(0.266,0.741,0.572),vec3(0.264,0.845,0.660)
);
// ----------------- Coastal:
// Jerlov 1C
const vec3 water_b_c_kd_1C[N_WLS] = vec3[](
    vec3(0.474,0.660,0.510),vec3(0.468,0.646,0.415),vec3(0.462,0.630,0.331),vec3(0.456,0.608,0.262),vec3(0.451,0.585,0.208),vec3(0.446,0.563,0.165),vec3(0.442,0.549,0.146),vec3(0.437,0.538,0.136),vec3(0.433,0.530,0.129),vec3(0.429,0.525,0.123),vec3(0.425,0.526,0.129),vec3(0.421,0.534,0.148),vec3(0.418,0.579,0.237),vec3(0.414,0.687,0.315),vec3(0.411,0.719,0.359),vec3(0.408,0.741,0.408),vec3(0.405,0.778,0.456),vec3(0.402,0.856,0.494),vec3(0.399,0.886,0.562),vec3(0.396,0.981,0.650)
);
// Jerlov 3C
const vec3 water_b_c_kd_3C[N_WLS] = vec3[](
    vec3(0.713,0.949,0.780),vec3(0.704,0.938,0.628),vec3(0.695,0.922,0.501),vec3(0.687,0.893,0.406),vec3(0.680,0.862,0.337),vec3(0.673,0.834,0.279),vec3(0.666,0.814,0.235),vec3(0.660,0.799,0.212),vec3(0.653,0.784,0.199),vec3(0.647,0.774,0.193),vec3(0.642,0.769,0.196),vec3(0.636,0.772,0.209),vec3(0.631,0.814,0.279),vec3(0.625,0.918,0.345),vec3(0.621,0.950,0.389),vec3(0.616,0.969,0.428),vec3(0.611,1.005,0.471),vec3(0.607,1.091,0.534),vec3(0.603,1.116,0.615),vec3(0.599,1.197,0.710)
);
// Jerlov 5C
const vec3 water_b_c_kd_5C[N_WLS] = vec3[](
    vec3(1.320,1.789,1.100),vec3(1.304,1.730,0.898),vec3(1.288,1.672,0.722),vec3(1.273,1.598,0.583),vec3(1.257,1.532,0.492),vec3(1.241,1.474,0.419),vec3(1.235,1.442,0.375),vec3(1.219,1.407,0.339),vec3(1.210,1.383,0.309),vec3(1.198,1.360,0.303),vec3(1.190,1.343,0.309),vec3(1.176,1.333,0.328),vec3(1.170,1.371,0.371),vec3(1.160,1.470,0.417),vec3(1.150,1.501,0.467),vec3(1.143,1.520,0.508),vec3(1.130,1.550,0.552),vec3(1.122,1.655,0.621),vec3(1.116,1.671,0.705),vec3(1.110,1.724,0.800));

// -----------------------------------------------------------------------------
// Returns the optical properties (scat b, ext c, diffuse ext kd) of the selected water at the given wl
/* Versions without interpolating waters, may be useful
// Evaluate a 1d->3d distribution (e.g., sensor response curve)
vec3 eval_3d(vec3 values[N_WLS], vec2 domain, float t) {
    float extent = domain.y-domain.x;
    float i_float = (t-domain.x)/extent; // 0 to 1
    i_float *= float(N_WLS); // 0.0 to N_WLS-1.0
    int i = int(floor(i_float)); // 0 to N_WLS-1
    // linear interp between the two neighboring ones:
    float i_fract = fract(i_float);
    return mix(values[i],values[i+1],i_fract); // i+1 may cause some trouble
}
// Evaluate a 1d->1d distribution (e.g. scattering spectrum)
float eval_1d(float values[N_WLS], vec2 domain, float t) {
    float extent = domain.y-domain.x;
    float i_float = (t-domain.x)/extent; // 0 to 1
    i_float *= float(N_WLS); // 0.0 to N_WLS-1.0
    int i = int(floor(i_float)); // 0 to N_WLS-1
    // linear interp between the two neighboring ones:
    float i_fract = fract(i_float);
    return mix(values[i],values[i+1],i_fract); // i+1 may cause some trouble
}
vec3 get_optical_props(float wl) {
    
    return eval_3d(water_b_c_kd, vec2(400.,700.), wl);
}
*/

// Type must be between 0. and 6. values in between interpolate the water tyoes
vec3 get_optical_props_anim(float wl, float type) {
    // Spectral sampling:
    vec2 domain = vec2(400.,700.);
    float extent = domain.y-domain.x;
    float i_float = (wl-domain.x)/extent; // 0 to 1
    i_float *= float(N_WLS); // 0.0 to N_WLS-1.0
    int i = int(floor(i_float)); // 0 to N_WLS-1
    // linear interp between the two neighboring ones:
    float i_fract = fract(i_float);
    
    // Type of water, 0 to 5:
    vec3 water1, water2;
    
    if (type<1.) {
        water1 = mix(water_b_c_kd_IB[i],water_b_c_kd_IB[i+1],i_fract);
        water2 = mix(water_b_c_kd_II[i],water_b_c_kd_II[i+1],i_fract);
    }
    else if (type<2.) {
        water1 = mix(water_b_c_kd_II[i],water_b_c_kd_II[i+1],i_fract);
        water2 = mix(water_b_c_kd_III[i],water_b_c_kd_III[i+1],i_fract);
    }
    else if (type<3.) {
        water1 = mix(water_b_c_kd_III[i],water_b_c_kd_III[i+1],i_fract);
        water2 = mix(water_b_c_kd_1C[i],water_b_c_kd_1C[i+1],i_fract);
    }
    else if (type<4.) {
        water1 = mix(water_b_c_kd_1C[i],water_b_c_kd_1C[i+1],i_fract);
        water2 = mix(water_b_c_kd_3C[i],water_b_c_kd_3C[i+1],i_fract);
    }
    else {
        water1 = mix(water_b_c_kd_3C[i],water_b_c_kd_3C[i+1],i_fract);
        water2 = mix(water_b_c_kd_5C[i],water_b_c_kd_5C[i+1],i_fract);
    }
    return mix(water1, water2, fract(type)); 
}

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// Shapes and scene from iq
// The MIT License
// Copyright © 2013 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// The license is here only not because I want to (can one
// license pieces of math?), but because people get upset
// if I don't add one...

// A list of useful distance function to simple primitives. All
// these functions (except for ellipsoid) return an exact
// euclidean distance, meaning they produce a better SDF than
// what you'd get if you were constructing them from boolean
// operations (such as cutting an infinite cylinder with two planes).

// List of other 3D SDFs:
//    https://www.shadertoy.com/playlist/43cXRl
// and
//    https://iquilezles.org/articles/distfunctions

//------------------------------------------------------------------
float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

float sdPlane( vec3 p )
{
	return p.y;
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sdBoxFrame( vec3 p, vec3 b, float e )
{
       p = abs(p  )-b;
  vec3 q = abs(p+e)-e;

  return min(min(
      length(max(vec3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
      length(max(vec3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
      length(max(vec3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}
float sdEllipsoid( in vec3 p, in vec3 r ) // approximated
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdTorus( vec3 p, vec2 t )
{
    return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

float sdCappedTorus(in vec3 p, in vec2 sc, in float ra, in float rb)
{
    p.x = abs(p.x);
    float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy,sc) : length(p.xy);
    return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}

float sdHexPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);

    const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec2 d = vec2(
       length(p.xy - vec2(clamp(p.x, -k.z*h.x, k.z*h.x), h.x))*sign(p.y - h.x),
       p.z-h.y );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdOctogonPrism( in vec3 p, in float r, float h )
{
  const vec3 k = vec3(-0.9238795325,   // sqrt(2+sqrt(2))/2 
                       0.3826834323,   // sqrt(2-sqrt(2))/2
                       0.4142135623 ); // sqrt(2)-1 
  // reflections
  p = abs(p);
  p.xy -= 2.0*min(dot(vec2( k.x,k.y),p.xy),0.0)*vec2( k.x,k.y);
  p.xy -= 2.0*min(dot(vec2(-k.x,k.y),p.xy),0.0)*vec2(-k.x,k.y);
  // polygon side
  p.xy -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  vec2 d = vec2( length(p.xy)*sign(p.y), p.z-h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}

float sdRoundCone( in vec3 p, in float r1, float r2, float h )
{
    vec2 q = vec2( length(p.xz), p.y );
    
    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(q,vec2(-b,a));
    
    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q-vec2(0.0,h)) - r2;
        
    return dot(q, vec2(a,b) ) - r1;
}

float sdRoundCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;
    
    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot2( pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float sdTriPrism( vec3 p, vec2 h )
{
    const float k = sqrt(3.0);
    h.x *= 0.5*k;
    p.xy /= h.x;
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x+k*p.y>0.0 ) p.xy=vec2(p.x-k*p.y,-k*p.x-p.y)/2.0;
    p.x -= clamp( p.x, -2.0, 0.0 );
    float d1 = length(p.xy)*sign(-p.y)*h.x;
    float d2 = abs(p.z)-h.y;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

// vertical
float sdCylinder( vec3 p, vec2 h )
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - h;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// arbitrary orientation
float sdCylinder(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);

    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}

// vertical
float sdCone( in vec3 p, in vec2 c, float h )
{
    vec2 q = h*vec2(c.x,-c.y)/c.y;
    vec2 w = vec2( length(p.xz), p.y );
    
	vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
    float k = sign( q.y );
    float d = min(dot( a, a ),dot(b, b));
    float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
	return sqrt(d)*sign(s);
}

float sdCappedCone( in vec3 p, in float h, in float r1, in float r2 )
{
    vec2 q = vec2( length(p.xz), p.y );
    
    vec2 k1 = vec2(r2,h);
    vec2 k2 = vec2(r2-r1,2.0*h);
    vec2 ca = vec2(q.x-min(q.x,(q.y < 0.0)?r1:r2), abs(q.y)-h);
    vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(dot2(ca),dot2(cb)) );
}

float sdCappedCone(vec3 p, vec3 a, vec3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;

    float x = sqrt( papa - paba*paba*baba );

    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;

    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );

    float cbx = x-ra - f*rba;
    float cby = paba - f;
    
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    
    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}

// c is the sin/cos of the desired cone angle
float sdSolidAngle(vec3 pos, vec2 c, float ra)
{
    vec2 p = vec2( length(pos.xz), pos.y );
    float l = length(p) - ra;
	float m = length(p - c*clamp(dot(p,c),0.0,ra) );
    return max(l,m*sign(c.y*p.x-c.x*p.y));
}

float sdOctahedron(vec3 p, float s)
{
    p = abs(p);
    float m = p.x + p.y + p.z - s;

    // exact distance
    #if 0
    vec3 o = min(3.0*p - m, 0.0);
    o = max(6.0*p - m*2.0 - o*3.0 + (o.x+o.y+o.z), 0.0);
    return length(p - s*o/(o.x+o.y+o.z));
    #endif
    
    // exact distance
    #if 1
 	vec3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;
    float k = clamp(0.5*(q.z-q.y+s),0.0,s); 
    return length(vec3(q.x,q.y-s+k,q.z-k)); 
    #endif
    
    // bound, not exact
    #if 0
	return m*0.57735027;
    #endif
}

float sdPyramid( in vec3 p, in float h )
{
    float m2 = h*h + 0.25;
    
    // symmetry
    p.xz = abs(p.xz);
    p.xz = (p.z>p.x) ? p.zx : p.xz;
    p.xz -= 0.5;
	
    // project into face plane (2D)
    vec3 q = vec3( p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y);
   
    float s = max(-q.x,0.0);
    float t = clamp( (q.y-0.5*p.z)/(m2+0.25), 0.0, 1.0 );
    
    float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
	float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);
    
    float d2 = min(q.y,-q.x*m2-q.y*0.5) > 0.0 ? 0.0 : min(a,b);
    
    // recover 3D and scale, and add sign
    return sqrt( (d2+q.z*q.z)/m2 ) * sign(max(q.z,-p.y));;
}

// la,lb=semi axis, h=height, ra=corner
float sdRhombus(vec3 p, float la, float lb, float h, float ra)
{
    p = abs(p);
    vec2 b = vec2(la,lb);
    float f = clamp( (ndot(b,b-2.0*p.xz))/dot(b,b), -1.0, 1.0 );
	vec2 q = vec2(length(p.xz-0.5*b*vec2(1.0-f,1.0+f))*sign(p.x*b.y+p.z*b.x-b.x*b.y)-ra, p.y-h);
    return min(max(q.x,q.y),0.0) + length(max(q,0.0));
}

float sdHorseshoe( in vec3 p, in vec2 c, in float r, in float le, vec2 w )
{
    p.x = abs(p.x);
    float l = length(p.xy);
    p.xy = mat2(-c.x, c.y, 
              c.y, c.x)*p.xy;
    p.xy = vec2((p.y>0.0 || p.x>0.0)?p.x:l*sign(-c.x),
                (p.x>0.0)?p.y:l );
    p.xy = vec2(p.x,abs(p.y-r))-vec2(le,0.0);
    
    vec2 q = vec2(length(max(p.xy,0.0)) + min(0.0,max(p.x,p.y)),p.z);
    vec2 d = abs(q) - w;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdU( in vec3 p, in float r, in float le, vec2 w )
{
    p.x = (p.y>0.0) ? abs(p.x) : length(p.xy);
    p.x = abs(p.x-r);
    p.y = p.y - le;
    float k = max(p.x,p.y);
    vec2 q = vec2( (k<0.0) ? -k : length(max(p.xy,0.0)), abs(p.z) ) - w;
    return length(max(q,0.0)) + min(max(q.x,q.y),0.0);
}


vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}


vec2 map( in vec3 pos )
{
    vec2 res = vec2( pos.y, 0.0 );

    // bounding box
    if( sdBox( pos-vec3(-2.0,0.3,0.25),vec3(0.3,0.3,1.0) )<res.x )
    {
      res = opU( res, vec2( sdSphere(    pos-vec3(-2.0,0.25, 0.0), 0.25 ), 26.9 ) );
	  res = opU( res, vec2( sdRhombus(  (pos-vec3(-2.0,0.25, 1.0)).xzy, 0.15, 0.25, 0.04, 0.08 ),17.0 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(0.0,0.3,-1.0),vec3(0.35,0.3,2.5) )<res.x )
    {
	res = opU( res, vec2( sdCappedTorus((pos-vec3( 0.0,0.30, 1.0))*vec3(1,-1,1), vec2(0.866025,-0.5), 0.25, 0.05), 25.0) );
    res = opU( res, vec2( sdBoxFrame(    pos-vec3( 0.0,0.25, 0.0), vec3(0.3,0.25,0.2), 0.025 ), 16.9 ) );
	res = opU( res, vec2( sdCone(        pos-vec3( 0.0,0.45,-1.0), vec2(0.6,0.8),0.45 ), 55.0 ) );
    res = opU( res, vec2( sdCappedCone(  pos-vec3( 0.0,0.25,-2.0), 0.25, 0.25, 0.1 ), 13.67 ) );
    res = opU( res, vec2( sdSolidAngle(  pos-vec3( 0.0,0.00,-3.0), vec2(3,4)/5.0, 0.4 ), 49.13 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(1.0,0.3,-1.0),vec3(0.35,0.3,2.5) )<res.x )
    {
	res = opU( res, vec2( sdTorus(      (pos-vec3( 1.0,0.30, 1.0)).xzy, vec2(0.25,0.05) ), 7.1 ) );
    res = opU( res, vec2( sdBox(         pos-vec3( 1.0,0.25, 0.0), vec3(0.3,0.25,0.1) ), 3.0 ) );
    res = opU( res, vec2( sdCapsule(     pos-vec3( 1.0,0.00,-1.0),vec3(-0.1,0.1,-0.1), vec3(0.2,0.4,0.2), 0.1  ), 31.9 ) );
	res = opU( res, vec2( sdCylinder(    pos-vec3( 1.0,0.25,-2.0), vec2(0.15,0.25) ), 8.0 ) );
    res = opU( res, vec2( sdHexPrism(    pos-vec3( 1.0,0.2,-3.0), vec2(0.2,0.05) ), 18.4 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(-1.0,0.35,-1.0),vec3(0.35,0.35,2.5))<res.x )
    {
	res = opU( res, vec2( sdPyramid(    pos-vec3(-1.0,-0.6,-3.0), 1.0 ), 13.56 ) );
	res = opU( res, vec2( sdOctahedron( pos-vec3(-1.0,0.15,-2.0), 0.35 ), 23.56 ) );
    res = opU( res, vec2( sdTriPrism(   pos-vec3(-1.0,0.15,-1.0), vec2(0.3,0.05) ),43.5 ) );
    res = opU( res, vec2( sdEllipsoid(  pos-vec3(-1.0,0.25, 0.0), vec3(0.2, 0.25, 0.05) ), 43.17 ) );
    res = opU( res, vec2( sdHorseshoe(  pos-vec3(-1.0,0.25, 1.0), vec2(cos(1.3),sin(1.3)), 0.2, 0.3, vec2(0.03,0.08) ), 11.5 ) );
    }

    // bounding box
    if( sdBox( pos-vec3(2.0,0.3,-1.0),vec3(0.35,0.3,2.5) )<res.x )
    {
    res = opU( res, vec2( sdOctogonPrism(pos-vec3( 2.0,0.2,-3.0), 0.2, 0.05), 51.8 ) );
    res = opU( res, vec2( sdCylinder(    pos-vec3( 2.0,0.14,-2.0), vec3(0.1,-0.1,0.0), vec3(-0.2,0.35,0.1), 0.08), 31.2 ) );
	res = opU( res, vec2( sdCappedCone(  pos-vec3( 2.0,0.09,-1.0), vec3(0.1,0.0,0.0), vec3(-0.2,0.40,0.1), 0.15, 0.05), 46.1 ) );
    res = opU( res, vec2( sdRoundCone(   pos-vec3( 2.0,0.15, 0.0), vec3(0.1,0.0,0.0), vec3(-0.1,0.35,0.1), 0.15, 0.05), 51.7 ) );
    res = opU( res, vec2( sdRoundCone(   pos-vec3( 2.0,0.20, 1.0), 0.2, 0.1, 0.3 ), 37.0 ) );
    }
    
    return res;
}

// https://iquilezles.org/articles/boxfunctions
vec2 iBox( in vec3 ro, in vec3 rd, in vec3 rad ) 
{
    vec3 m = 1.0/rd;
    vec3 n = m*ro;
    vec3 k = abs(m)*rad;
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;
	return vec2( max( max( t1.x, t1.y ), t1.z ),
	             min( min( t2.x, t2.y ), t2.z ) );
}

vec2 raycast( in vec3 ro, in vec3 rd )
{
    vec2 res = vec2(-1.0,-1.0);

    float tmin = 1.0;
    float tmax = 20.0;
    // raytrace floor plane
    #ifndef NO_FLOOR
    float tp1 = (0.0-ro.y)/rd.y;
    if( tp1>0.0 )
    {
        tmax = min( tmax, tp1 );
        res = vec2( tp1, 1.0 );
    }
    #endif
    
    // raymarch primitives   
    vec2 tb = iBox( ro-vec3(0.0,0.4,-0.5), rd, vec3(2.5,0.41,3.0) );
    if( tb.x<tb.y && tb.y>0.0 && tb.x<tmax)
    {
        //return vec2(tb.x,2.0);
        tmin = max(tb.x,tmin);
        tmax = min(tb.y,tmax);

        float t = tmin;
        for( int i=0; i<70 && t<tmax; i++ )
        {
            vec2 h = map( ro+rd*t );
            if( abs(h.x)<(0.0001*t) )
            { 
                res = vec2(t,h.y); 
                break;
            }
            t += h.x;
        }
    }
    
    return res;
}

// https://iquilezles.org/articles/rmshadows
float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
    // bounding volume
    float tp = (0.8-ro.y)/rd.y; if( tp>0.0 ) tmax = min( tmax, tp );

    float res = 1.0;
    float t = mint;
    for( int i=0; i<24; i++ )
    {
		float h = map( ro + rd*t ).x;
        float s = clamp(8.0*h/t,0.0,1.0);
        res = min( res, s );
        t += clamp( h, 0.01, 0.2 );
        if( res<0.004 || t>tmax ) break;
    }
    res = clamp( res, 0.0, 1.0 );
    return res*res*(3.0-2.0*res);
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos )
{
#if 0
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*map( pos + e.xyy ).x + 
					  e.yyx*map( pos + e.yyx ).x + 
					  e.yxy*map( pos + e.yxy ).x + 
					  e.xxx*map( pos + e.xxx ).x );
#else
    // inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
    vec3 n = vec3(0.0);
    for( int i=0; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+0.0005*e).x;
      //if( n.x+n.y+n.z>100.0 ) break;
    }
    return normalize(n);
#endif    
}

// https://iquilezles.org/articles/nvscene2008/rwwtt.pdf
float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float h = 0.01 + 0.12*float(i)/4.0;
        float d = map( pos + h*nor ).x;
        occ += (h-d)*sca;
        sca *= 0.95;
        if( occ>0.35 ) break;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 ) * (0.5+0.5*nor.y);
}

// https://iquilezles.org/articles/checkerfiltering
float checkersGradBox( in vec2 p, in vec2 dpdx, in vec2 dpdy )
{
    // filter kernel
    vec2 w = abs(dpdx)+abs(dpdy) + 0.001;
    // analytical integral (box filter)
    vec2 i = 2.0*(abs(fract((p-0.5*w)*0.5)-0.5)-abs(fract((p+0.5*w)*0.5)-0.5))/w;
    // xor pattern
    return 0.5 - 0.5*i.x*i.y;                  
}

mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
{
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.0);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv =          ( cross(cu,cw) );
    return mat3( cu, cv, cw );
}


// End of iq's code

//------------------------------------------------------------------
//------------------------------------------------------------------

// Spectra to colors:


// -------------------- CIE XYZ 2006 ----------------------------

// Source http://cvrl.ioo.ucl.ac.uk/index.htm:

const vec2 CIEXYZ_SPAN = vec2(380.0, 830.0);
const int CIEXYZ_SAMPLES = 89;
const vec3 ciexyz[CIEXYZ_SAMPLES] = vec3[](
    vec3(3.769647E-03,4.146161E-04,1.847260E-02),
    vec3(9.382967E-03,1.059646E-03,4.609784E-02),
    vec3(2.214302E-02,2.452194E-03,1.096090E-01),
    vec3(4.742986E-02,4.971717E-03,2.369246E-01),
    vec3(8.953803E-02,9.079860E-03,4.508369E-01),
    vec3(1.446214E-01,1.429377E-02,7.378822E-01),
    vec3(2.035729E-01,2.027369E-02,1.051821E+00),
    vec3(2.488523E-01,2.612106E-02,1.305008E+00),
    vec3(2.918246E-01,3.319038E-02,1.552826E+00),
    vec3(3.227087E-01,4.157940E-02,1.748280E+00),
    vec3(3.482554E-01,5.033657E-02,1.917479E+00),
    vec3(3.418483E-01,5.743393E-02,1.918437E+00),
    vec3(3.224637E-01,6.472352E-02,1.848545E+00),
    vec3(2.826646E-01,7.238339E-02,1.664439E+00),
    vec3(2.485254E-01,8.514816E-02,1.522157E+00),
    vec3(2.219781E-01,1.060145E-01,1.428440E+00),
    vec3(1.806905E-01,1.298957E-01,1.250610E+00),
    vec3(1.291920E-01,1.535066E-01,9.991789E-01),
    vec3(8.182895E-02,1.788048E-01,7.552379E-01),
    vec3(4.600865E-02,2.064828E-01,5.617313E-01),
    vec3(2.083981E-02,2.379160E-01,4.099313E-01),
    vec3(7.097731E-03,2.850680E-01,3.105939E-01),
    vec3(2.461588E-03,3.483536E-01,2.376753E-01),
    vec3(3.649178E-03,4.277595E-01,1.720018E-01),
    vec3(1.556989E-02,5.204972E-01,1.176796E-01),
    vec3(4.315171E-02,6.206256E-01,8.283548E-02),
    vec3(7.962917E-02,7.180890E-01,5.650407E-02),
    vec3(1.268468E-01,7.946448E-01,3.751912E-02),
    vec3(1.818026E-01,8.575799E-01,2.438164E-02),
    vec3(2.405015E-01,9.071347E-01,1.566174E-02),
    vec3(3.098117E-01,9.544675E-01,9.846470E-03),
    vec3(3.804244E-01,9.814106E-01,6.131421E-03),
    vec3(4.494206E-01,9.890228E-01,3.790291E-03),
    vec3(5.280233E-01,9.994608E-01,2.327186E-03),
    vec3(6.133784E-01,9.967737E-01,1.432128E-03),
    vec3(7.016774E-01,9.902549E-01,8.822531E-04),
    vec3(7.967750E-01,9.732611E-01,5.452416E-04),
    vec3(8.853376E-01,9.424569E-01,3.386739E-04),
    vec3(9.638388E-01,8.963613E-01,2.117772E-04),
    vec3(1.051011E+00,8.587203E-01,1.335031E-04),
    vec3(1.109767E+00,8.115868E-01,8.494468E-05),
    vec3(1.143620E+00,7.544785E-01,5.460706E-05),
    vec3(1.151033E+00,6.918553E-01,3.549661E-05),
    vec3(1.134757E+00,6.270066E-01,2.334738E-05),
    vec3(1.083928E+00,5.583746E-01,1.554631E-05),
    vec3(1.007344E+00,4.895950E-01,1.048387E-05),
    vec3(9.142877E-01,4.229897E-01,0.000000E+00),
    vec3(8.135565E-01,3.609245E-01,0.000000E+00),
    vec3(6.924717E-01,2.980865E-01,0.000000E+00),
    vec3(5.755410E-01,2.416902E-01,0.000000E+00),
    vec3(4.731224E-01,1.943124E-01,0.000000E+00),
    vec3(3.844986E-01,1.547397E-01,0.000000E+00),
    vec3(2.997374E-01,1.193120E-01,0.000000E+00),
    vec3(2.277792E-01,8.979594E-02,0.000000E+00),
    vec3(1.707914E-01,6.671045E-02,0.000000E+00),
    vec3(1.263808E-01,4.899699E-02,0.000000E+00),
    vec3(9.224597E-02,3.559982E-02,0.000000E+00),
    vec3(6.639960E-02,2.554223E-02,0.000000E+00),
    vec3(4.710606E-02,1.807939E-02,0.000000E+00),
    vec3(3.292138E-02,1.261573E-02,0.000000E+00),
    vec3(2.262306E-02,8.661284E-03,0.000000E+00),
    vec3(1.575417E-02,6.027677E-03,0.000000E+00),
    vec3(1.096778E-02,4.195941E-03,0.000000E+00),
    vec3(7.608750E-03,2.910864E-03,0.000000E+00),
    vec3(5.214608E-03,1.995557E-03,0.000000E+00),
    vec3(3.569452E-03,1.367022E-03,0.000000E+00),
    vec3(2.464821E-03,9.447269E-04,0.000000E+00),
    vec3(1.703876E-03,6.537050E-04,0.000000E+00),
    vec3(1.186238E-03,4.555970E-04,0.000000E+00),
    vec3(8.269535E-04,3.179738E-04,0.000000E+00),
    vec3(5.758303E-04,2.217445E-04,0.000000E+00),
    vec3(4.058303E-04,1.565566E-04,0.000000E+00),
    vec3(2.856577E-04,1.103928E-04,0.000000E+00),
    vec3(2.021853E-04,7.827442E-05,0.000000E+00),
    vec3(1.438270E-04,5.578862E-05,0.000000E+00),
    vec3(1.024685E-04,3.981884E-05,0.000000E+00),
    vec3(7.347551E-05,2.860175E-05,0.000000E+00),
    vec3(5.259870E-05,2.051259E-05,0.000000E+00),
    vec3(3.806114E-05,1.487243E-05,0.000000E+00),
    vec3(2.758222E-05,1.080001E-05,0.000000E+00),
    vec3(2.004122E-05,7.863920E-06,0.000000E+00),
    vec3(1.458792E-05,5.736935E-06,0.000000E+00),
    vec3(1.068141E-05,4.211597E-06,0.000000E+00),
    vec3(7.857521E-06,3.106561E-06,0.000000E+00),
    vec3(5.768284E-06,2.286786E-06,0.000000E+00),
    vec3(4.259166E-06,1.693147E-06,0.000000E+00),
    vec3(3.167765E-06,1.262556E-06,0.000000E+00),
    vec3(2.358723E-06,9.422514E-07,0.000000E+00),
    vec3(1.762465E-06,7.053860E-07,0.000000E+00)
);


vec3 xyz2srgb(vec3 xyz) {
    // From https://github.com/tobspr/GLSL-Color-Spaces
    const mat3 XYZ_2_RGB = (mat3(
         3.2404542,-0.9692660, 0.0556434,
        -1.5371385, 1.8760108,-0.2040259,
        -0.4985314, 0.0415560, 1.0572252
    ));
    return XYZ_2_RGB * xyz;
}


// aces (from https://64.github.io/tonemapping/#aces)
vec3 aces_approx(vec3 v)
{
    //v *= 0.6f;
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((v*(a*v+b))/(v*(c*v+d)+e), 0.0f, 1.0f);
}

// converts a linear srgb value to gamma (the pipeline is usually 
// render->inner_product(spectra,cmf)->XYZ->linear_sRGB->(tonemap here->linear_sRGB_LDR)->gamma_sRGB )
vec3 gamma_correct(vec3 linear_srgb)
{
    vec3 a = 12.92 * linear_srgb;
    vec3 b = 1.055 * pow(linear_srgb, vec3(1.0 / 2.4)) - 0.055;
    vec3 c = step(vec3(0.0031308), linear_srgb);
    return mix(a, b, c);
}



// Returns the XYZ sensitivity at the given wl
// (note that this is likely quite silly, I don't think we would have the same
//  sensitivity opening our eyes underwater)
vec3 sensitivity(float wl) {
    float extent = CIEXYZ_SPAN.y-CIEXYZ_SPAN.x;
    float i_float = (wl-CIEXYZ_SPAN.x)/extent; // 0 to 1
    i_float *= float(CIEXYZ_SAMPLES); // 0.0 to N_WLS-1.0
    int i = int(floor(i_float)); // 0 to N_WLS-1
    // linear interp between the two neighboring ones:
    float i_fract = fract(i_float);
    return mix(ciexyz[i],ciexyz[i+1],i_fract); 
}
// ^ end of COMMON
// -----------------------------------------------------------



// -----------------------------------------------------------
// MAIN

// This shader is an example application implementing Real-Time Underwater Spectral Rendering 
// [https://doi.org/10.1111/cgf.15009]
// The scene geometry and base lighting is a fork of iq's https://www.shadertoy.com/view/Xds3zN

// The method renders a number of natural waters using measured spectral coefficients (see Common tab)
// (scattering, extinction and diffuse downwelling attenuation coefficient, Kd)
// You can set the Jerlov water type at the start of this file (0 to 5, from least to most turbid)
// by default, it changes with time

// Here the scene geometry is raymarched, but the inscattering only needs the scene distance, so 
// it can be used in a rasterization pipeline

// For spatially varying volumetric single scattering, in the paper we used froxel lighting.
// Here, we use raymarching instead (slower but simpler). 
// ------------------------------------------------------------------------------
// NOTE: I disabled SS by default but I recommend to enable it if you can (line 43)!
// ------------------------------------------------------------------------------

// More credits:
// The caustics here are from Dave Hoskins (https://www.shadertoy.com/view/MdlXz8)
// The dithering (very important to avoid banding) is from https://www.shadertoy.com/view/4dcSRX and 
// https://www.shadertoy.com/view/Md3SRX, which were adapted from valve's GDC15 talk by alex vlachos

//------------------------------------------------------------------
// The water types are 0.0 to 5.0 (three open, three coastal, from less turbid to more turbid)
// Uncomment to animate the water type:
#define ANIMATE_WATER_TYPE
const float JERLOV_WATER_IDX = 0.; // if the above is commented, this sets the water type (0 to 5)

// Alternatively the depth of the scene changes w/time
// #define ANIMATE_DEPTH 

//------------------------------------------------------------------
// Performance:

// Volumetric SS:
// 0-> Disabled (cheapest..)
// 1-> Brute force
#if HW_PERFORMANCE==0
#define SINGLE_SCATTERING 0
#else
// --------------------------------------------------------------
#define SINGLE_SCATTERING 0    // change this to 1 if you can
// --------------------------------------------------------------
#endif
// Single scattering (volumetric shadows&caustics etc): 
#if SINGLE_SCATTERING == 1
// You can tune the SS settings here:
// A ray will be marched from the origin to the minimum of intersection point and ss_dist, 
// evaluating lighting at ss_steps steps along the way (it would be nice to do this at a quarter res or so)
const int ss_steps = 40; // can be increased to reduce some artifacts... for a price
const float ss_dist = 4.; // maximum range for ss
const float ss_multiplier = 0.0004; // multiplier for the strength of ss
const float ss_caustics_mult = 1.; // multiplier for the caustics in ss
#endif

//------------------------------------------------------------------
// 0-> no caustics
// 1-> one caustic texture
// 2-> two caustic textures (nicer, more cost) (animated differently to hide repetitions)
// 3-> 2 and, if SS is enabled, it also shades the volumetric caustics
#if HW_PERFORMANCE==0
#define CAUSTICS 1 
#else
#define CAUSTICS 3   
#endif

//------------------------------------------------------------------
// Scene lighting settings:

// scene depth (m below the surface):
#ifndef ANIMATE_DEPTH
#define SCENE_DEPTH 25. 
#else
#define SCENE_DEPTH 25. + 5. + 10.*sin(iTime* 1.)
//5.*sin(iTime)
#endif
const float turbidity = .9; // Multiplier of the water coefs (10x coefs is the same as the whole scene being 10x scale)

const vec3 sun_dir = vec3(-0.5, 1.4, -0.6);
const float sun_power = 40.;
const float exposure = 2.5;

// Spectral sampling (fast diminishing results):
const int samples = 8;

const float MAX_DIST_INSCATTERING = 1000.; // only for the analytical approx - no cost 

#if CAUSTICS > 0
const float caustic_power = 9.0;
const int caustic_iters = 4;
const float caustic_speed = 1.1;
const float caustic_freq = 1.2; // spatial freq of caustics
#endif


#define SPLIT_CHECK if (uv.x>iMouse.x/iResolution.x) 


//------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Implements Eq 9 https://doi.org/10.1111/cgf.15009
// (Note: unlike in the paper here y points up)
// samples: nb of spectral samples
// range: min and max wls
// irradiance_in: the irradiance at the water surface
// o: origin of the ray (water surface assumed to be at the origin)
// w: direction of the ray, must be normalized
// dist: distance to the intersection (inscattering will be computed between o and o+w*dist)
// (Note: it would be nice to multiply the spectra with D65, the daylight spd. Right now it is 
//  illuminated by a constant spectrum (float irradiance_in). This works as a cloudy day I guess)
vec3 integrate_L_anim(int samples, vec2 range, float irradiance_in,
                vec3 o, vec3 w, float dist, float turbidity,
                float water_type) {
    float o_depth = -o.y; // depth of the origin underwater
    float cosine = -w.y;  // cosine of the direction (positive pointing down) 
    
    vec3 col = vec3(0.);
    for (int i=0;i<samples;i++) {
        float x = (float(i)+0.5)/float(samples); // important to center the samples if using few
        float wl = mix(range.x, range.y, x);
        
        vec3 props = turbidity * get_optical_props_anim(wl, water_type); // scat, ext, kd
        float b = props.x, c = props.y, kd = props.z; // Unpack the coefs
        
        float floor_refl = 0.; // TODO: add reflectance from floor (now infinitely deep.. or completely black)
        float c_mod = c + kd*cosine; // modified c with kd:
        
        // Eq 9 in the paper (sigma_s=b; sigma_t=c):
        float L = (b * irradiance_in) * (1. - exp(-dist*c_mod)) 
                  * (exp(-kd*o_depth) + floor_refl) / (4.*PI*c_mod);
                  
        vec3 sens = sensitivity(wl);
        col += L*sens;
    }
    col /= float(samples);
    return col;
}

const float TAU = 6.28318530718;

// Very nice tileable caustics
// from @DaveHoskins https://www.shadertoy.com/view/MdlXz8
float caustic(in vec2 p, in float time) {
	time = time * .5+23.0;

	vec2 i = vec2(p);
	float c = 1.0;
	float inten = .005;
    
	for (int n = 0; n < caustic_iters; n++) 
	{
		float t = time * (1.0 - (3.5 / float(n+1)));
		i = p + vec2(cos(t - i.x) + sin(t + i.y), sin(t - i.y) + cos(t + i.x));
		c += 1.0/length(vec2(p.x / (sin(i.x+t)/inten),p.y / (cos(i.y+t)/inten)));
	}
	c /= float(caustic_iters);
	c = 1.17-pow(c, 1.4);
    return pow(abs(c), 8.0);
}


float getCaustics(vec3 p, vec3 sunDir) {
    float t_surf = (10.0-p.y)/sunDir.y; // dist to water surface
    if( t_surf<=0.0 ) return 0.; // above water?
    vec3 p_surf = p + sunDir*t_surf;
    p_surf *= 0.11 * caustic_freq;
    
    vec2 uv =  mod(iTime+p_surf.xz*TAU*2.0, TAU)-250.0; 
    float c = caustic(uv, iTime*caustic_speed);     
#if CAUSTICS > 1
    // add a second texture, offset and with slightly diff parameters: 
    uv =  mod(-vec2(0.4,0.7)*iTime+p_surf.xz*TAU*2.0*0.56, TAU)-250.0;  
    c += 1.1 * caustic(uv, iTime*caustic_speed*1.1); 
#endif
    return c * caustic_power;
}


// This just returns an id between 0. and 5. which is one of the water types.
// (Float values in between interpolate the two neighboring waters, e.g., 0.5 is a mix of Jerlov IB and II)
// If ANIMATE_WATER_TYPE is set, it loops between them
// otherwise this returns JERLOV_WATER_IDX
float get_water_type() {
    #ifndef ANIMATE_WATER_TYPE
    return JERLOV_WATER_IDX;
    #else
    // note: removed the last water here bc it is very very turbid
    //return 5.-abs(mod((iTime*0.6), 10.)-5.); 
    return 4.-abs((mod((iTime*0.6), 8.)-4.)); 
    #endif
}

// -----------------------------------------------------------------

// AA from iq's shader:
#define ZERO (min(iFrame,0))
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 1   // make this 2 or 3 for antialiasing
#endif


// -------------------------------------------------------------------------------------------------
// Dithering:
// Dithering is important because water produces very smooth gradients that are quite visible with
// 255 bits. You can turn it off (or try DITHER_TEST) to check it out (espetially visible in more turbid waters)
#define DITHER 1
// 0-> none; 1-> https://www.shadertoy.com/view/Md3SRX

// For comparison, enable this so that mouse.x controls where dither is enabled:
//#define DITHER_TEST 

// from valve's GDC15 talk by alex vlachos, implemented by 
// https://www.shadertoy.com/view/4dcSRX and https://www.shadertoy.com/view/Md3SRX
vec3 ScreenSpaceDither(vec2 vScreenPos)
{
    vec3 vDither = vec3(dot(vec2(131.0, 312.0), vScreenPos.xy + iTime));
    vDither.rgb = fract(vDither.rgb / vec3(103.0, 71.0, 97.0)) - vec3(0.5, 0.5, 0.5);
    return 0.05*(vDither.rgb / 15.0);
}


vec2 intersectWaterSurface(vec3 ro, vec3 ld) {
    vec2 res = vec2(-1.0,-1.0); // (distance, id={-1 not found, 1 found water )
    float tp1 = (0.0-ro.y)/-ld.y;
    if( tp1>0.0 )
    {
        res = vec2( tp1, 1.0 );
    }
    return res;
}

float upsample(vec3 col, float wl) {
    //Low budget version of spectral upsampling (something like this [https://rgl.epfl.ch/publications/Jakob2019Spectral] would be better)
    if (wl<480.) return col.b; 
    else if (wl<560.) return col.g;
    else return col.r;
}

// Calculate the full SS at ro in direction rd, up to t units away
// for a light from ld
vec3 calcSS(in vec3 ro, in vec3 rd, in vec3 ld, float t, float turbidity, vec2 uv) {
    #if SINGLE_SCATTERING == 1
    t = min(t,ss_dist);
    float water_type = get_water_type();
    vec3 l_ss = vec3(0.);
    // ss is already expensive, so dont make it spectral:
    float ext = get_optical_props_anim(550., water_type).y;
    
    for (int i=1; i<ss_steps+1; i++) {
        float t_step = (-0.5 + t * float(i)) / float(ss_steps);
        
        // vec3 p = mix(ro, ro+rd*t, float(i)/float(ss_steps));
        vec3 pos = ro+ rd*t_step;

        float shadow = calcSoftshadow( pos, ld, 0.02, 2.5 );
        vec2 water_surf = intersectWaterSurface( pos, ld );
        if (water_surf.y<0.) water_surf.y *=-1.;//return 0.*vec3(1.,0.,1.); // Not found
        
        float t_surface = water_surf.x; // Distance to water surface
        
        float l_caustic = 1.;
        #if CAUSTICS > 2
        // get caustics here for nice godrays
        l_caustic = 0.2 * getCaustics(pos, ld);
        l_caustic = ss_caustics_mult*1.*pow(l_caustic, 5.);
        #endif
        
        if (shadow>0.2) {
            float transmittance = exp(-(t_surface + t_step) * ext); // extinction
            l_ss += vec3(0.4,0.6,0.8) * shadow * (sun_power*l_caustic) 
                  * transmittance;         
        }
    }
    return l_ss * ss_multiplier /float(ss_steps);
    
    #endif
    return vec3(0.,0.,0.);
}

vec3 render( in vec3 ro, in vec3 rd, in vec3 rdx, in vec3 rdy, vec2 uv )
{ 

    float water_type = get_water_type();
    // make the last waters less turbid for visualization in this scene:
    float turbidity_mod = mix(turbidity, turbidity/6., water_type/6.);
   
    // background
    vec3 col = vec3(0.);
        
    // raycast scene
    vec2 res = raycast(ro,rd);
    float t = res.x; // dist
	float m = res.y; // Material
    
    if (m>-0.5) {
    	vec3 lin = vec3(0.0);

        vec3 pos = ro + t*rd;
        vec3 nor = (m<1.5) ? vec3(0.0,1.0,0.0) : calcNormal( pos );
        vec3 ref = reflect( rd, nor );
        
        // material        
        col = 0.2 + 0.2*sin( m*2.0 + vec3(0.0,1.0,2.0) );
        float ks = 1.0;
        
        if( m<1.5 )
        {
            // project pixel footprint into the plane
            vec3 dpdx = ro.y*(rd/rd.y-rdx/rdx.y);
            vec3 dpdy = ro.y*(rd/rd.y-rdy/rdy.y);

            float f = checkersGradBox( 3.0*pos.xz, 3.0*dpdx.xz, 3.0*dpdy.xz );
            col = 0.15 + f*vec3(0.05);
            ks = 0.4;
        }

        // lighting
        float occ = calcAO( pos, nor );
        
        { // surface shading, adapted from iq's
            float dif = sqrt(clamp( 0.5+0.5*nor.y, 0.0, 1.0 ));
                  dif *= occ;
            float spe = smoothstep( -0.2, 0.2, ref.y );
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0+dot(nor,rd),0.0,1.0), 5.0 );
                  spe *= calcSoftshadow( pos, ref, 0.02, 2.5 );
            #if CAUSTICS > 0
            vec3  lig = normalize( sun_dir );
            //vec3  hal = normalize( lig-rd );
            float dif2 = clamp( dot( nor, lig ), 0.0, 1.0 );
                  dif2 *= calcSoftshadow( pos, lig, 0.02, 2.5 );
            float c = getCaustics(pos, lig);
            dif += dif2 * c * 0.8;
            spe = spe * (1. + c*c*0.4);//sharpened specular caustics
            #endif
            
            float pos_depth = -pos.y + SCENE_DEPTH;

            float wl_step = (700.-400.)/float(samples);
            for (float wl=400.+wl_step*0.5; wl<700.; wl+=wl_step) {
                vec3 b_c_kd = turbidity_mod*get_optical_props_anim(wl, water_type); // turbidity is a multiplier for artistic control
                float l_wl = (upsample(col, wl) * 0.9 * dif  // Diffuse
                             + .9*spe*ks)
                             * exp(-b_c_kd.z * pos_depth) * // diffuse extinction (kd) from water surface to pos
                               exp(-b_c_kd.y * t );      // extinction (c) from pos to ro (t units)
                lin += sensitivity(wl) * l_wl;
            }
            lin /= float(samples);
        }
       col = lin;
    }

    // -------------------------------------------------------------
    // water
    t = min(t, MAX_DIST_INSCATTERING);
    
    // MS lighting:
    col += integrate_L_anim(samples, vec2(400.,700.), sun_power*1., 
            ro-vec3(0.,SCENE_DEPTH,0.), 
            normalize(rd), t, 
            turbidity_mod, water_type);
            
    // Single scattering, if enabled:
    col += calcSS(ro, normalize(rd), normalize(sun_dir), t, turbidity_mod, uv );
    
	return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 mo = iMouse.xy/iResolution.xy;
	float time = 32.0 + iTime*1.5;

    // camera	
    
    // move target
    vec3 ta = vec3( 0.25, mix(1.6, -1.05, 1.-mo.y), -0.75 );
    vec3 ro = ta + vec3( 4.5*cos(0.1*time + 7.0*mo.x), 2.4, 4.5*sin(0.1*time + 7.0*mo.x) );
    ro.y = mix(ro.y, ta.y+2.2, mo.y);
    // camera-to-world transformation
    mat3 ca = setCamera( ro, ta, 0. );

    vec3 tot = vec3(0.0);
#if AA>1
    for( int m=ZERO; m<AA; m++ )
    for( int n=ZERO; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (2.0*(fragCoord+o)-iResolution.xy)/iResolution.y;
#else    
        vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
#endif
        // focal length
        const float fl = 2.5;
        
        // ray direction
        vec3 rd = ca * normalize( vec3(p,fl) );

         // ray differentials
        vec2 px = (2.0*(fragCoord+vec2(1.0,0.0))-iResolution.xy)/iResolution.y;
        vec2 py = (2.0*(fragCoord+vec2(0.0,1.0))-iResolution.xy)/iResolution.y;
        vec3 rdx = ca * normalize( vec3(px,fl) );
        vec3 rdy = ca * normalize( vec3(py,fl) );
        
        // render	
        vec3 col = render( ro, rd, rdx, rdy, fragCoord.xy/iResolution.xy);

        // gain
        col = col*3.0/(2.5+col);
        
		// gamma
        //col =pow( col, vec3(0.4545) );
        //col = gamma_correct(exposure*col); // apply gain for this one
        col = aces_approx(exposure*col);

        tot += col;
        
#if AA>1
    }
    tot /= float(AA*AA);
#endif

#if DITHER > 0
    // Dither to avoid banding (see https://www.shadertoy.com/view/XdtXWn)
    // could be done to the params instead (e.g. inscattering distance or angle)
    #ifdef DITHER_TEST
    vec2 uv = fragCoord.xy/iResolution.xy;
    bool dither = uv.x < iMouse.x/iResolution.x;
    tot *= smoothstep(0.,1.,abs(uv.x-iMouse.x/iResolution.x)*1000.);
    #else
    bool dither = true;
    #endif
    if (dither) tot += ScreenSpaceDither(fragCoord);
#endif
    
    fragColor = vec4( tot, 1.0 );
}
