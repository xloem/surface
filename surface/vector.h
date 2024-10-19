#include <tgmath.h>

// probably want at some point a transformation structure
// maybe a 4x3 matrix alternatively a quat and a pt3
// likely redo the functions below named 'transform' at that point

typedef float real;

typedef union {
    struct { real x, y; };
    struct { real u, v; };
    struct { real r, i; };
    real idx[2];
} vec2;
typedef vec2 pt2;
typedef vec2 cplx;
typedef vec2 cplx2;

typedef union {
    struct { real x, y, z; };
    vec2 xy;
    real idx[3];
} vec3;
typedef vec3 pt3;

typedef union {
    struct { real w; vec3 xyz; };
    struct { real r; union { vec3 ijk; struct{ real i, j, k; }; }; };
    struct { real a, b, c, d; };
    real idx[4];
} vec4;
typedef vec4 quat;
typedef vec4 vers;
typedef vec4 cplx4;

#define AXES2 (vec2[2]){\
    {1,0},              \
    {0,1}               \
}
#define NEGAXES2 (vec2[2]){\
    {-1,0},                \
    {0,-1}                 \
}
#define U AXES2[0]
#define V AXES2[1]
#define R2 AXES2[0]
#define I2 AXES2[1]
#define NEGU NEGAXES2[0]
#define NEGV NEGAXES2[1]
#define NEGR2 NEGAXES2[0]
#define NEGI2 NEGAXES2[1]

#define AXES3 (vec3[3]){\
    {1,0,0},            \
    {0,1,0},            \
    {0,0,1}             \
}
#define NEGAXES3 (vec3[3]){\
    {-1,0,0},              \
    {0,-1,0},              \
    {0,0,-1}               \
}
#define X AXES3[0]
#define Y AXES3[1]
#define Z AXES3[2]
#define NEGX NEGAXES3[0]
#define NEGY NEGAXES3[1]
#define NEGZ NEGAXES3[2]

#define AXES4 (vec4[4]){\
    {1,0,0,0},          \
    {0,1,0,0},          \
    {0,0,1,0},          \
    {0,0,0,1}           \
}
#define NEGAXES4 (vec4[4]){\
    {-1,0,0,0},            \
    {0,-1,0,0},            \
    {0,0,-1,0},            \
    {0,0,0,-1}             \
}
#define R4 AXES4[0]
#define I4 AXES4[1]
#define J4 AXES4[2]
#define K4 AXES4[3]
#define NEGR4 NEGAXES4[0]
#define NEGI4 NEGAXES4[1]
#define NEGJ4 NEGAXES4[2]
#define NEGK4 NEGAXES4[3]

#define _def_kernel_ov(name, o_t, odecl, v_t, vdecl, kernel) \
extern inline o_t* name##into(o_t* _o, v_t _v)               \
{                                                            \
    const int odims = sizeof(o_t)/sizeof(real);              \
    const int vdims = sizeof(v_t)/sizeof(real);              \
    const int dims = vdims > odims ? vdims : odims;          \
    odecl; vdecl;                                            \
    for (int _dim = 0; _dim < dims; ++ _dim)                 \
    {                                                        \
        int odim = odims > 1 ? _dim : 0;                     \
        int vdim = vdims > 1 ? _dim : 0;                     \
        kernel;                                              \
    }                                                        \
    return _o;                                               \
}                                                            \
extern inline o_t name(v_t _v)                               \
{ o_t _o; return *name##into(&_o, _v); }                     \

#define _def_kernel_olr(name, o_t, odecl, l_t, ldecl, r_t, rdecl, kernel)\
extern inline o_t* name##into(o_t* _o, l_t _l, r_t _r)                   \
{                                                                        \
    const int odims = sizeof(o_t)/sizeof(real);                          \
    const int ldims = sizeof(l_t)/sizeof(real);                          \
    const int rdims = sizeof(r_t)/sizeof(real);                          \
    const int dims = ldims > rdims ? ldims : rdims;                      \
    odecl; ldecl; rdecl;                                                 \
    for (int _dim = 0; _dim < dims; ++ _dim)                             \
    {                                                                    \
        int odim = odims > 1 ? _dim : 0;                                 \
        int ldim = ldims > 1 ? _dim : 0;                                 \
        int rdim = rdims > 1 ? _dim : 0;                                 \
        kernel;                                                          \
    }                                                                    \
    return _o;                                                           \
}                                                                        \
extern inline o_t name(l_t _l, r_t _r)                                   \
{ o_t _o; return *name##into(&_o, _l, _r); }                             \

#define def_kernel_out1valN(name, init, kernel) \
_def_kernel_ov(name##2,                         \
    real,  real* out = _o; *out = init,         \
    vec2,  real* val = _v.idx,                  \
    kernel)                                     \
_def_kernel_ov(name##3,                         \
    real,  real* out = _o; *out = init,         \
    vec3,  real* val = _v.idx,                  \
    kernel)                                     \
_def_kernel_ov(name##4,                         \
    real,  real* out = _o; *out = init,         \
    vec4,  real* val = _v.idx,                  \
    kernel)                                     \

#define def_kernel_outNvalN(name, kernel)  \
_def_kernel_ov(name##2,                    \
    vec2,  real* out = _o->idx,            \
    vec2,  real* val = _v.idx,             \
    kernel)                                \
_def_kernel_ov(name##3,                    \
    vec3,  real* out = _o->idx,            \
    vec3,  real* val = _v.idx,             \
    kernel)                                \
_def_kernel_ov(name##4,                    \
    vec4,  real* out = _o->idx,            \
    vec4,  real* val = _v.idx,             \
    kernel)                                \

#define def_kernel_out1leftNrightN(name, init, kernel)\
_def_kernel_olr(name##2,                              \
    real,  real* out = _o; *out = init,               \
    vec2,  real* left = _l.idx,                       \
    vec2,  real* right = _r.idx,                      \
    kernel)                                           \
_def_kernel_olr(name##3,                              \
    real,  real* out = _o; *out = init,               \
    vec3,  real* left = _l.idx,                       \
    vec3,  real* right = _r.idx,                      \
    kernel)                                           \
_def_kernel_olr(name##4,                              \
    real,  real* out = _o; *out = init,               \
    vec3,  real* left = _l.idx,                       \
    vec3,  real* right = _r.idx,                      \
    kernel)                                           \

#define def_kernel_outNleftNrightN(name, kernel)\
_def_kernel_olr(                                \
    name##2,                                    \
    vec2,  real* out = _o->idx,                 \
    vec2,  real* left = _l.idx,                 \
    vec2,  real* right = _r.idx,                \
    kernel)                                     \
_def_kernel_olr(                                \
    name##21,                                   \
    vec2,  real* out = _o->idx,                 \
    vec2,  real* left = _l.idx,                 \
    real,  real* right = &_r,                   \
    kernel)                                     \
_def_kernel_olr(                                \
    name##12,                                   \
    vec2,  real* out = _o->idx,                 \
    real,  real* left = &_l,                    \
    vec2,  real* right = _r.idx,                \
    kernel)                                     \
_def_kernel_olr(                                \
    name##3,                                    \
    vec3,  real* out = _o->idx,                 \
    vec3,  real* left = _l.idx,                 \
    vec3,  real* right = _r.idx,                \
    kernel)                                     \
_def_kernel_olr(                                \
    name##31,                                   \
    vec3,  real* out = _o->idx,                 \
    vec3,  real* left = _l.idx,                 \
    real,  real* right = &_r,                   \
    kernel)                                     \
_def_kernel_olr(                                \
    name##13,                                   \
    vec3,  real* out = _o->idx,                 \
    real,  real* left = &_l,                    \
    vec3,  real* right = _r.idx,                \
    kernel)                                     \
_def_kernel_olr(                                \
    name##4,                                    \
    vec4,  real* out = _o->idx,                 \
    vec4,  real* left = _l.idx,                 \
    vec4,  real* right = _r.idx,                \
    kernel)                                     \
_def_kernel_olr(                                \
    name##41,                                   \
    vec4,  real* out = _o->idx,                 \
    vec4,  real* left = _l.idx,                 \
    real,  real* right = &_r,                   \
    kernel)                                     \
_def_kernel_olr(                                \
    name##14,                                   \
    vec4,  real* out = _o->idx,                 \
    real,  real* left = &_l,                    \
    vec4,  real* right = _r.idx,                \
    kernel)                                     \

def_kernel_out1valN(normsq, 0, *out += val[vdim] * val[vdim]);
#define norm2(val2) sqrt(normsq2(val2))
#define norm3(val3) sqrt(normsq3(val3))
#define norm4(val4) sqrt(normsq4(val4))

def_kernel_outNvalN(neg, out[odim] = -val[vdim])
def_kernel_outNvalN(conj, out[odim] = vdim ? -val[vdim] : val[vdim])
def_kernel_outNvalN(negconj, out[odim] = vdim ? val[vdim] : -val[vdim])

def_kernel_out1leftNrightN(dot, 0, *out += left[ldim] * right[rdim])

def_kernel_outNleftNrightN(add, out[odim] = left[ldim] + right[rdim])
def_kernel_outNleftNrightN(sub, out[odim] = left[ldim] - right[rdim])
def_kernel_outNleftNrightN(mul, out[odim] = left[ldim] * right[rdim])
def_kernel_outNleftNrightN(div, out[odim] = left[ldim] / right[rdim])

#define normed2into(out, val2) div2into(out, val2, norm2(val2))
#define normed3into(out, val3) div3into(out, val3, norm3(val3))
#define normed4into(out, val4) div4into(out, val4, norm4(val4))
#define normed2(val2) div2(val2, norm2(val2))
#define normed3(val3) div3(val3, norm3(val3))
#define normed4(val4) div4(val4, norm4(val4);

extern inline cplx* mulcplx2into(cplx* out, cplx left, cplx right)
{
    out->r = left.r*right.r - left.i*right.i;
    out->i = left.r*right.i + left.i*right.r;
}
extern inline cplx mulcplx2(cplx left, cplx right)
{ cplx out; return *mulcplx2into(&out, left, right); }

extern inline pt3* crossinto(pt3* out, pt3 left, pt3 right)
{
    out->x = left.y*right.z - left.z*right.y;
    out->y = left.z*right.x - left.x*right.z;
    out->z = left.x*right.y - left.y*right.x;
    return out;
}
extern inline pt3 cross(pt3 left, pt3 right)
{ pt3 out; return *crossinto(&out, left, right); }

extern inline quat* mulcplx4into(quat* out, quat left, quat right)
{
    out->r= left.r*right.r - left.i*right.i - left.j*right.j - left.k*right.k;
    out->i= left.r*right.i + left.i*right.r + left.j*right.k - left.k*right.j;
    out->j= left.r*right.j - left.i*right.k + left.j*right.r + left.k*right.i;
    out->k= left.r*right.k + left.i*right.j - left.j*right.i + left.k*right.r;
    return out;
}
extern inline quat mulcplx4(quat left, quat right)
{ quat out; return *mulcplx4into(&out, left, right); }

extern inline quat* inv4into(quat* out, quat val)
{
    return div41into(out, *conj4into(out, val), normsq4(val));
}
extern inline quat inv4(quat val)
{ quat out; return *inv4into(&out, val); }
extern inline vers* invversinto(vers* out, vers val)
{
    return conj4into(out, val);
}
extern inline vers invvers(vers val)
{ vers out; return *invversinto(&out, val); }

extern inline quat* logversinto(quat* out, vers val)
{
    out->r = 0;
    mul31into(&out->ijk, val.ijk, acos(val.r)/norm3(val.ijk));
    return out;
}
extern inline quat logvers(vers val)
{ quat out; return *logversinto(&out, val); }
extern inline quat* logquatinto(quat* out, quat val)
{
    real normsq = normsq4(val);
    if (normsq != 1) {
        real norm = sqrt(normsq);
        out->r = log(norm);
        mul31into(&out->ijk, val.ijk, acos(val.r/norm)/norm3(val.ijk));
        return out;
    } else {
        return logversinto(out, val);
    }
}
extern inline quat logquat(quat val)
{ quat out; return *logquatinto(&out, val); }
extern inline quat* expquatinto(quat* out, quat val)
{
    real norm = norm3(val.ijk);
    out->r = cos(norm);
    mul31into(&out->ijk, val.ijk, sin(norm)/norm);
    return mul41into(out, *out, exp(val.r));
}
extern inline quat expquat(quat val)
{ quat out; return *expquatinto(&out, val); }
#define expversinto expquatinto
#define expvers expquat
extern inline quat* powquatinto(quat* out, quat val, real exp)
{
    return expquatinto(out, *mul41into(out, *logquatinto(out, val), exp));
}
extern inline quat powquat(quat val, real exp)
{ quat out; return *powquatinto(&out, val, exp); }
extern inline vers* powversinto(vers* out, quat val, real exp)
{
    return expversinto(out, *mul41into(out, *logversinto(out, val), exp));
}
extern inline vers powvers(vers val, real exp)
{ vers out; return *powversinto(&out, val, exp); }

extern inline vers* fromaxisangleinto(vers* out, real radians, pt3 axis)
{
    real halfradians = radians / 2;
    real axis_normsq = normsq3(axis);
    if (axis_normsq != 1) {
        div31into(&axis, axis, sqrt(axis_normsq));
    }
    out->r = cos(halfradians);
    mul31into(&out->ijk, axis, sin(halfradians));
    return out;
}
extern inline vers fromaxisangle(real radians, pt3 axis)
{ vers out; return *fromaxisangleinto(&out, radians, axis); }

#define negversinto negconj4into
#define negvers negconj4
#define negquatinto negconj4into
#define negquat negconj4

#define mulversinto mulcplx4into
#define mulvers mulcplx4
#define mulquatinto mulcplx4into
#define mulquat mulcplx4

#define rotate4into(out, first, second) mulversinto(out, second, first)
#define rotate4(first, second) mulvers(second, first)

extern inline pt3* rotateinto(pt3* out, pt3 pt, vers orientation)
{
    // rot.ijk = sin(theta/2) * u
    // rot.r = cos(theta/2)
    // cross product distributes with multiplication
    // so sin(theta/2) * (u cross p) = rot.ijk cross p
    // r = p + 2 * rot.r * cross3(rot.ijk, p) + 2 * cross3(rot.ijk, cross3(rot.ijk, p))
    // out = vec + 2 * rot.r * cross + 2 * cross3(rot.ijk, cross)
    vec3 cross;
    crossinto(&cross, orientation.ijk, pt);
    crossinto(out, orientation.ijk, cross);
    add3into(out, *out, mul31(cross, orientation.r));
    mul31into(out, *out, 2);
    add3into(out, *out, pt);
    return out;
}
extern inline pt3 rotate(pt3 vec, vers orientation)
{ pt3 out; return *rotateinto(&out, vec, orientation); }
extern inline pt3 transform(pt3 vec, vers orientation)
{
    vec4 p = { .r = 0, .ijk = vec };
    mul4into(&p, orientation, p);
    mul4into(&p, p, inv4(orientation));
    return p.xyz;
}
extern inline pt3* transforminto(pt3* out, pt3 vec, vers orientation)
{ *out = transform(vec, orientation); return out; }

extern inline pt3(*transformationaxesinto(pt3(*out)[3], quat orientation))[3]
{
    real s = 1 / normsq4(orientation);
    real bs, cs, ds;
    if (s != 1) {
        s = 2 * sqrt(s);
        bs = orientation.b*s; cs = orientation.c*s; ds = orientation.d*s;
    } else {
        bs = orientation.b*2; cs = orientation.c*2; ds = orientation.d*2;
    }
    real ab = orientation.a*bs, ac = orientation.a*cs, ad = orientation.a*ds;
    real bb = orientation.b*bs, bc = orientation.b*cs, bd = orientation.b*ds;
    real cc = orientation.c*cs, cd = orientation.c*ds, dd = orientation.d*ds;
    (*out)[0] = (pt3){1-cc-dd,  bc-ad,  bd+ac};
    (*out)[1] = (pt3){  bc+ad,1-bb-dd,  cd-ab};
    (*out)[2] = (pt3){  bd-ac,  cd+ab,1-bb-cc};
    return out;
}

extern inline pt3(*rotationaxesinto(pt3(*out)[3], vers orientation))[3]
{
    real s = 1 / normsq4(orientation);
    real bs, cs, ds;
    if (s != 1) {
        s = 2 * s;
        bs = orientation.b*s; cs = orientation.c*s; ds = orientation.d*s;
    } else {
        bs = orientation.b*2; cs = orientation.c*2; ds = orientation.d*2;
    }
    real ab = orientation.a*bs, ac = orientation.a*cs, ad = orientation.a*ds;
    real bb = orientation.b*bs, bc = orientation.b*cs, bd = orientation.b*ds;
    real cc = orientation.c*cs, cd = orientation.c*ds, dd = orientation.d*ds;
    (*out)[0] = (pt3){1-cc-dd,  bc-ad,  bd+ac};
    (*out)[1] = (pt3){  bc+ad,1-bb-dd,  cd-ab};
    (*out)[2] = (pt3){  bd-ac,  cd+ab,1-bb-cc};
    return out;
}
/* to convert a precise matrix to a versor see https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion */

extern inline vec4* toaxisangleinto(vec4* out, vers orientation)
{
    real norm = norm3(orientation.ijk);
    div31into(&out->xyz, orientation.ijk, norm);
    out->w = 2 * atan2(norm, orientation.r);
    return out;
}
extern inline vec4 toaxisangle(vers orientation)
{ vec4 out; return *toaxisangleinto(&out, orientation); }

extern inline vers* slerpinto(vers* out, vers left, vers right, real t)
{
    // slerp 
    // = qo * pow(inv(q0)*q1, t)
    // = q1 * pow(inv(q1)*q0, 1-t)
    // = pow(q0 * inv(q1), 1-t) * q1
    // = pow(q1 * inv(q0), t) * q0
    mulversinto(out, *invversinto(out, left), right);
    return mulversinto(out, left, *powversinto(out, *out, t));
}
extern inline vers slerp(vers left, vers right, real t)
{ vers out; return *slerpinto(&out, left, right, t); }

extern inline vers* slerpdifferentialinto(vers* out, vers left, vers right, real t)
{
    return logversinto(out, *mulversinto(out, right, *invversinto(out, left)));
}
extern inline vers slerpdifferential(vers left, vers right, real t)
{ vers out; return *slerpdifferentialinto(&out, left, right, t); }

extern inline vec4(*rotationdifferentialinto(vec4(*out)[4], vec3 pt, vers orientation))[4]
{
    quat pq = mulquat((quat){.w = 0, .xyz = pt}, orientation);
    sub4into(out[0], pq, conj4(pq));
    quat pqi = mulquat(pq, I4);
    sub4into(out[1], conj4(pqi), pqi);
    quat pqj = mulquat(pq, J4);
    sub4into(out[2], conj4(pqj), pqj);
    quat pqk = mulquat(pq, K4);
    sub4into(out[3], conj4(pqk), pqk);
    return out;
}
