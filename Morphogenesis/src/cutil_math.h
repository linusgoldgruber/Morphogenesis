    /*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/*
    This file implements common mathematical operations on vector types
    (cl_float3, float4 etc.) since these are not provided as standard by CUDA.

    The syntax is modelled on the Cg standard library.
*/

//#ifndef CUTIL_MATH_H
//#define CUTIL_MATH_H
//
#ifndef CL_UTIL_MATH_H
#define CL_UTIL_MATH_H

//#include "cuda_runtime.h"
//
#include "CL/cl.h"

////////////////////////////////////////////////////////////////////////////////
typedef unsigned int uint;
typedef unsigned short ushort;

/*
#ifndef __CUDACC__
#include <math.h>

inline float fminf(float a, float b)
{
  return a < b ? a : b;
}

inline float fmaxf(float a, float b)
{
  return a > b ? a : b;
}

inline int max(int a, int b)
{
  return a > b ? a : b;
}

inline int min(int a, int b)
{
  return a < b ? a : b;
}
#endif
*/

#ifndef __OPENCLC__
#include <math.h>

inline float fminf(float a, float b)
{
  return a < b ? a : b;
}

inline float fmaxf(float a, float b)
{
  return a > b ? a : b;
}

inline int max(int a, int b)
{
  return a > b ? a : b;
}

inline int min(int a, int b)
{
  return a < b ? a : b;
}
#endif

// float functions
////////////////////////////////////////////////////////////////////////////////

// lerp
inline __device__ __host__ float lerp(float a, float b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float clamp(float f, float a, float b)
{
    return fmaxf(a, fminf(f, b));
}

// int2 functions
////////////////////////////////////////////////////////////////////////////////

// negate
inline __host__ __device__ int2 operator-(int2 &a)
{
    return make_int2(-a.x, -a.y);
}

// addition
inline __host__ __device__ int2 operator+(int2 a, int2 b)
{
    return make_int2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(int2 &a, int2 b)
{
    a.x += b.x; a.y += b.y;
}

// subtract
inline __host__ __device__ int2 operator-(int2 a, int2 b)
{
    return make_int2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(int2 &a, int2 b)
{
    a.x -= b.x; a.y -= b.y;
}

// multiply
inline __host__ __device__ int2 operator*(int2 a, int2 b)
{
    return make_int2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ int2 operator*(int2 a, int s)
{
    return make_int2(a.x * s, a.y * s);
}
inline __host__ __device__ int2 operator*(int s, int2 a)
{
    return make_int2(a.x * s, a.y * s);
}
inline __host__ __device__ void operator*=(int2 &a, int s)
{
    a.x *= s; a.y *= s;
}

// float2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float2 make_float2(float s)
{
    return make_float2(s, s);
}
inline __host__ __device__ float2 make_float2(int2 a)
{
    return make_float2(float(a.x), float(a.y));
}

// negate
inline __host__ __device__ float2 operator-(float2 &a)
{
    return make_float2(-a.x, -a.y);
}

// addition
inline __host__ __device__ float2 operator+(float2 a, float2 b)
{
    return make_float2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(float2 &a, float2 b)
{
    a.x += b.x; a.y += b.y;
}

// subtract
inline __host__ __device__ float2 operator-(float2 a, float2 b)
{
    return make_float2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(float2 &a, float2 b)
{
    a.x -= b.x; a.y -= b.y;
}

// multiply
inline __host__ __device__ float2 operator*(float2 a, float2 b)
{
    return make_float2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ float2 operator*(float2 a, float s)
{
    return make_float2(a.x * s, a.y * s);
}
inline __host__ __device__ float2 operator*(float s, float2 a)
{
    return make_float2(a.x * s, a.y * s);
}
inline __host__ __device__ void operator*=(float2 &a, float s)
{
    a.x *= s; a.y *= s;
}

// divide
inline __host__ __device__ float2 operator/(float2 a, float2 b)
{
    return make_float2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ float2 operator/(float2 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ float2 operator/(float s, float2 a)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ void operator/=(float2 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// lerp
inline __device__ __host__ float2 lerp(float2 a, float2 b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float2 clamp(float2 v, float a, float b)
{
    return make_float2(clamp(v.x, a, b), clamp(v.y, a, b));
}

inline __device__ __host__ float2 clamp(float2 v, float2 a, float2 b)
{
    return make_float2(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y));
}

// dot product
inline __host__ __device__ float dot(float2 a, float2 b)
{ 
    return a.x * b.x + a.y * b.y;
}

// length
inline __host__ __device__ float length(float2 v)
{
    return sqrtf(dot(v, v));
}

// normalize
inline __host__ __device__ float2 normalize(float2 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float2 floor(const float2 v)
{
    return make_float2(floor(v.x), floor(v.y));
}

// reflect
inline __host__ __device__ float2 reflect(float2 i, float2 n)
{
	return i - 2.0f * n * dot(n,i);
}

// cl_float3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ cl_float3 make_cl_float3(float s)
{
    return make_cl_float3(s, s, s);
}
inline __host__ __device__ cl_float3 make_cl_float3(float2 a)
{
    return make_cl_float3(a.x, a.y, 0.0f);
}
inline __host__ __device__ cl_float3 make_cl_float3(float2 a, float s)
{
    return make_cl_float3(a.x, a.y, s);
}
inline __host__ __device__ cl_float3 make_cl_float3(float4 a)
{
    return make_cl_float3(a.x, a.y, a.z);  // discards w
}
inline __host__ __device__ cl_float3 make_cl_float3(cl_int3 a)
{
    return make_cl_float3(float(a.x), float(a.y), float(a.z));
}

// negate
inline __host__ __device__ cl_float3 operator-(cl_float3 &a)
{
    return make_cl_float3(-a.x, -a.y, -a.z);
}

// min
static __inline__ __host__ __device__ cl_float3 fminf(cl_float3 a, cl_float3 b)
{
	return make_cl_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}

// max
static __inline__ __host__ __device__ cl_float3 fmaxf(cl_float3 a, cl_float3 b)
{
	return make_cl_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}

// addition
inline __host__ __device__ cl_float3 operator+(cl_float3 a, cl_float3 b)
{
    return make_cl_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ cl_float3 operator+(cl_float3 a, float b)
{
    return make_cl_float3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(cl_float3 &a, cl_float3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ cl_float3 operator-(cl_float3 a, cl_float3 b)
{
    return make_cl_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ cl_float3 operator-(cl_float3 a, float b)
{
    return make_cl_float3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(cl_float3 &a, cl_float3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ cl_float3 operator*(cl_float3 a, cl_float3 b)
{
    return make_cl_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ cl_float3 operator*(cl_float3 a, float s)
{
    return make_cl_float3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ cl_float3 operator*(float s, cl_float3 a)
{
    return make_cl_float3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(cl_float3 &a, float s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ cl_float3 operator/(cl_float3 a, cl_float3 b)
{
    return make_cl_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ cl_float3 operator/(cl_float3 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ cl_float3 operator/(float s, cl_float3 a)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ void operator/=(cl_float3 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// lerp
inline __device__ __host__ cl_float3 lerp(cl_float3 a, cl_float3 b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ cl_float3 clamp(cl_float3 v, float a, float b)
{
    return make_cl_float3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ cl_float3 clamp(cl_float3 v, cl_float3 a, cl_float3 b)
{
    return make_cl_float3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}

// dot product
inline __host__ __device__ float dot(cl_float3 a, cl_float3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// cross product
inline __host__ __device__ cl_float3 cross(cl_float3 a, cl_float3 b)
{ 
    return make_cl_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); 
}

// length
inline __host__ __device__ float length(cl_float3 v)
{
    return sqrtf(dot(v, v));
}

// normalize
inline __host__ __device__ cl_float3 normalize(cl_float3 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ cl_float3 floor(const cl_float3 v)
{
    return make_cl_float3(floor(v.x), floor(v.y), floor(v.z));
}

// reflect
inline __host__ __device__ cl_float3 reflect(cl_float3 i, cl_float3 n)
{
	return i - 2.0f * n * dot(n,i);
}

// float4 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float4 make_float4(float s)
{
    return make_float4(s, s, s, s);
}
inline __host__ __device__ float4 make_float4(cl_float3 a)
{
    return make_float4(a.x, a.y, a.z, 0.0f);
}
inline __host__ __device__ float4 make_float4(cl_float3 a, float w)
{
    return make_float4(a.x, a.y, a.z, w);
}
inline __host__ __device__ float4 make_float4(int4 a)
{
    return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));
}

// negate
inline __host__ __device__ float4 operator-(float4 &a)
{
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}

// min
static __inline__ __host__ __device__ float4 fminf(float4 a, float4 b)
{
	return make_float4(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z), fminf(a.w,b.w));
}

// max
static __inline__ __host__ __device__ float4 fmaxf(float4 a, float4 b)
{
	return make_float4(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z), fmaxf(a.w,b.w));
}

// addition
inline __host__ __device__ float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w);
}
inline __host__ __device__ void operator+=(float4 &a, float4 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

// subtract
inline __host__ __device__ float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w);
}
inline __host__ __device__ void operator-=(float4 &a, float4 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

// multiply
inline __host__ __device__ float4 operator*(float4 a, float s)
{
    return make_float4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ float4 operator*(float s, float4 a)
{
    return make_float4(a.x * s, a.y * s, a.z * s, a.w * s);
}
inline __host__ __device__ void operator*=(float4 &a, float s)
{
    a.x *= s; a.y *= s; a.z *= s; a.w *= s;
}

// divide
inline __host__ __device__ float4 operator/(float4 a, float4 b)
{
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}
inline __host__ __device__ float4 operator/(float4 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ float4 operator/(float s, float4 a)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ void operator/=(float4 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// lerp
inline __device__ __host__ float4 lerp(float4 a, float4 b, float t)
{
    return a + t*(b-a);
}

// clamp
inline __device__ __host__ float4 clamp(float4 v, float a, float b)
{
    return make_float4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b));
}

inline __device__ __host__ float4 clamp(float4 v, float4 a, float4 b)
{
    return make_float4(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z), clamp(v.w, a.w, b.w));
}

// dot product
inline __host__ __device__ float dot(float4 a, float4 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// length
inline __host__ __device__ float length(float4 r)
{
    return sqrtf(dot(r, r));
}

// normalize
inline __host__ __device__ float4 normalize(float4 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float4 floor(const float4 v)
{
    return make_float4(floor(v.x), floor(v.y), floor(v.z), floor(v.w));
}

// cl_int3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ cl_int3 make_cl_int3(int s)
{
    return make_cl_int3(s, s, s);
}
inline __host__ __device__ cl_int3 make_cl_int3(cl_float3 a)
{
    return make_cl_int3(int(a.x), int(a.y), int(a.z));
}

// negate
inline __host__ __device__ cl_int3 operator-(cl_int3 &a)
{
    return make_cl_int3(-a.x, -a.y, -a.z);
}

// min
inline __host__ __device__ cl_int3 min(cl_int3 a, cl_int3 b)
{
    return make_cl_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}

// max
inline __host__ __device__ cl_int3 max(cl_int3 a, cl_int3 b)
{
    return make_cl_int3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}

// addition
inline __host__ __device__ cl_int3 operator+(cl_int3 a, cl_int3 b)
{
    return make_cl_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(cl_int3 &a, cl_int3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ cl_int3 operator-(cl_int3 a, cl_int3 b)
{
    return make_cl_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ void operator-=(cl_int3 &a, cl_int3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ cl_int3 operator*(cl_int3 a, cl_int3 b)
{
    return make_cl_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ cl_int3 operator*(cl_int3 a, int s)
{
    return make_cl_int3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ cl_int3 operator*(int s, cl_int3 a)
{
    return make_cl_int3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(cl_int3 &a, int s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ cl_int3 operator/(cl_int3 a, cl_int3 b)
{
    return make_cl_int3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ cl_int3 operator/(cl_int3 a, int s)
{
    return make_cl_int3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ cl_int3 operator/(int s, cl_int3 a)
{
    return make_cl_int3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ void operator/=(cl_int3 &a, int s)
{
    a.x /= s; a.y /= s; a.z /= s;
}

// clamp
inline __device__ __host__ int clamp(int f, int a, int b)
{
    return max(a, min(f, b));
}

inline __device__ __host__ cl_int3 clamp(cl_int3 v, int a, int b)
{
    return make_cl_int3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ cl_int3 clamp(cl_int3 v, cl_int3 a, cl_int3 b)
{
    return make_cl_int3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}


// ucl_int3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ ucl_int3 make_ucl_int3(uint s)
{
    return make_ucl_int3(s, s, s);
}
inline __host__ __device__ ucl_int3 make_ucl_int3(cl_float3 a)
{
    return make_ucl_int3(uint(a.x), uint(a.y), uint(a.z));
}

// min
inline __host__ __device__ ucl_int3 min(ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(min(a.x,b.x), min(a.y,b.y), min(a.z,b.z));
}

// max
inline __host__ __device__ ucl_int3 max(ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(max(a.x,b.x), max(a.y,b.y), max(a.z,b.z));
}

// addition
inline __host__ __device__ ucl_int3 operator+(ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ void operator+=(ucl_int3 &a, ucl_int3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ ucl_int3 operator-(ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ void operator-=(ucl_int3 &a, ucl_int3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ ucl_int3 operator*(ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ ucl_int3 operator*(ucl_int3 a, uint s)
{
    return make_ucl_int3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ ucl_int3 operator*(uint s, ucl_int3 a)
{
    return make_ucl_int3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(ucl_int3 &a, uint s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ ucl_int3 operator/(ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ ucl_int3 operator/(ucl_int3 a, uint s)
{
    return make_ucl_int3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ ucl_int3 operator/(uint s, ucl_int3 a)
{
    return make_ucl_int3(a.x / s, a.y / s, a.z / s);
}
inline __host__ __device__ void operator/=(ucl_int3 &a, uint s)
{
    a.x /= s; a.y /= s; a.z /= s;
}

// clamp
inline __device__ __host__ uint clamp(uint f, uint a, uint b)
{
    return max(a, min(f, b));
}

inline __device__ __host__ ucl_int3 clamp(ucl_int3 v, uint a, uint b)
{
    return make_ucl_int3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ ucl_int3 clamp(ucl_int3 v, ucl_int3 a, ucl_int3 b)
{
    return make_ucl_int3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}

#endif
