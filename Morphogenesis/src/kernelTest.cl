//----------------------------------------------------------------------------------
//
// FLUIDS v.3 - SPH Fluid Simulator for CPU and GPU
// Copyright (C) 2012-2013. Rama Hoetzlein, http://fluids3.com
//
// BSD 3-clause:
// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this 
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this 
//    list of conditions and the following disclaimer in the documentation and/or 
//    other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may 
//    be used to endorse or promote products derived from this software without specific 
//   prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//----------------------------------------------------------------------------------

//#include </home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/fluid_system_opencl.cl>
//#include <float.h>        //No
//#include <stdint.h>       //No
//#include <string.h>       //No
#include <assert.h>         //Yes
#include </home/goldi/Documents/Libraries/RandomCL/generators/well512.cl>       //Yes
//#include "/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/fluid.h"        //No
#include <CL/cl.h>        //No



// __constant FParams		fparam;			// CPU Fluid params
// __constant FBufs			fbuf;			// GPU Particle buffers (unsorted). An FBufs struct holds an array of pointers.
// __constant FBufs			ftemp;			// GPU Particle buffers (sorted)
// __constant FGenome		fgenome;		// GPU Genome for particle automata behaviour. Also holds morphogen diffusability.
// __constant_FBondParams    fbondparams;    // GPU copy of remodelling parameters.
// __constant uint			gridActive;
//
#define SCAN_BLOCKSIZE		512
//#define FLT_MIN  0.000000001              // set here as 2^(-30)
//#define UINT_MAX 65535
//
// //if(fparam.debug>2) => device printf
//
//
// __kernel void insertParticlesCLTest ( int pnum ) {
//     uint i = get_global_id(0);
//     if ( i >= pnum ) return;
//
//     float3  gridMin     =  fparam.gridMin;
//     float3  gridDelta   =  fparam.gridDelta;
//     int3    gridRes     =  fparam.gridRes;
//     int3    gridScan    =  fparam.gridScanMax;
//     int     gridTot     =  fparam.gridTotal;
//
//     int     gs;
//     float3  gcf;
//     int3    gc;
//
//     gcf = (fbuf.bufF3(FPOS)[i] - gridMin) * gridDelta;
//     gc  = (int3)(gcf.x, gcf.y, gcf.z);
//     gs  = (gc.y * gridRes.z + gc.z) * gridRes.x + gc.x;
//
//     if ( gc.x >= 1 && gc.x <= gridScan.x && gc.y >= 1 && gc.y <= gridScan.y && gc.z >= 1 && gc.z <= gridScan.z ) {
//         fbuf.bufI(FGCELL)[i] = gs;
//         fbuf.bufI(FGNDX)[i] = atomic_add(&fbuf.bufI(FGRIDCNT)[gs], 1);
//
//         for(int gene=0; gene<NUM_GENES; gene++){
//             if (fbuf.bufI(FEPIGEN)[i + gene*fparam.maxPoints] > 0){
//                 atomic_add(&fbuf.bufI(FGRIDCNT_ACTIVE_GENES)[gene*gridTot + gs], 1);
//             }
//         }
//     } else {
//         fbuf.bufI(FGCELL)[i] = GRID_UNDEF;
//     }
// }

__kernel void prefixFixupTest(__global uint *input, __global uint *aux, int len) {

    unsigned int t = get_local_id(0);
    unsigned int start = t + 2 * get_group_id(0) * SCAN_BLOCKSIZE;
    if (start < len)					input[start] += aux[get_group_id(0)];
    if (start + SCAN_BLOCKSIZE < len)   input[start + SCAN_BLOCKSIZE] += aux[get_group_id(0)];
}
