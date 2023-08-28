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

#define CL_KERNEL
#define SCAN_BLOCKSIZE		512
#define GRID_UNDEF			4294967295
//#define FLT_MIN  0.000000001              // set here as 2^(-30) //REDEFINED
//#define UINT_MAX 65535    //REDEFINED

#include "fluid_system_opencl.h"
#include "fluid.h"
#include </home/goldi/Documents/Libraries/RandomCL/generators/well512.cl>


__constant struct FParams	fparam  = {
        .debug = 0,
		.numThreads = 0, .numBlocks = 0, .threadsPerBlock = 0,
		.gridThreads = 0, .gridBlocks = 0,
		.szPnts = 0,
		.szGrid = 0,
		.stride = 0, .pnum = 0, .pnumActive = 0, .maxPoints = 0,
        .freeze = 0,
        .frame = 0,
		.chk = 0,
		.pdist = 0, .pmass = 0, .prest_dens = 0,
		.pextstiff = 0, .pintstiff = 0,
		.pradius = 0, .psmoothradius = 0, .r2 = 0, .psimscale = 0, .pvisc = 0, .psurface_t = 0,
		.pforce_min = 0, .pforce_max = 0, .pforce_freq = 0, .pground_slope = 0,
		.pvel_limit = 0, .paccel_limit = 0, .pdamp = 0,
		.pboundmin = 0, .pboundmax = 0, .pgravity = 0,
		.AL = 0, .AL2 = 0, .VL = 0, .VL2 = 0,
		.H = 0, .d2 = 0, .rd2 = 0, .vterm = 0,		// used in force calculation
		.poly6kern = 0, .spikykern = 0, .lapkern = 0, .gausskern = 0, .wendlandC2kern = 0,
		.gridSize = 0, .gridDelta = 0, .gridMin = 0, .gridMax = 0,
		.gridRes = 0, .gridScanMax = 0,
		.gridSrch = 0, .gridTotal = 0, .gridAdjCnt = 0, .gridActive = 0,
        .actuation_factor = 0,
        .actuation_period = 0,
		.gridAdj = {0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0}
	};			// CPU Fluid params
__constant struct FBufs			fbuf = {};
			// GPU Particle buffers (unsorted). An FBufs struct holds an array of pointers.
__constant struct FBufs			ftemp = {};			// GPU Particle buffers (sorted)

__constant struct FGenome		fgenome = {};		// GPU Genome for particle automata behaviour. Also holds morphogen diffusability.
//__constant_FBondParams    fbondparams;    // GPU copy of remodelling parameters.
//__constant uint			gridActive;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// Kernels /////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__kernel void prefixFixup(__global uint *input, __global uint *aux, int len) {

    unsigned int t = get_local_id(0);
    unsigned int start = t + 2 * get_group_id(0) * SCAN_BLOCKSIZE;
    if (start < len)					input[start] += aux[get_group_id(0)];
    if (start + SCAN_BLOCKSIZE < len)   input[start + SCAN_BLOCKSIZE] += aux[get_group_id(0)];
}

__kernel void prefixSum(__global uint *input, __global uint *output, __global uint *aux, int len, int zeroff)
{
    __local uint scan_array[SCAN_BLOCKSIZE << 1];
    unsigned int t1 = get_local_id(0) + 2 * get_group_id(0) * SCAN_BLOCKSIZE;
    unsigned int t2 = t1 + SCAN_BLOCKSIZE;

    // Pre-load into shared memory
    scan_array[get_local_id(0)] = (t1<len) ? input[t1] : 0.0f;
    scan_array[get_local_id(0) + SCAN_BLOCKSIZE] = (t2<len) ? input[t2] : 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Reduction
    int stride;
    for (stride = 1; stride <= SCAN_BLOCKSIZE; stride <<= 1) {
        int index = (get_local_id(0) + 1) * stride * 2 - 1;
        if (index < 2 * SCAN_BLOCKSIZE)
            scan_array[index] += scan_array[index - stride];
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Post reduction
    for (stride = SCAN_BLOCKSIZE >> 1; stride > 0; stride >>= 1) {
        int index = (get_local_id(0) + 1) * stride * 2 - 1;
        if (index + stride < 2 * SCAN_BLOCKSIZE)
            scan_array[index + stride] += scan_array[index];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Output values & aux
    if (t1 + zeroff < len) output[t1 + zeroff] = scan_array[get_local_id(0)];
    if (t2 + zeroff < len) output[t2 + zeroff] = (get_local_id(0) == SCAN_BLOCKSIZE - 1 && zeroff) ? 0 : scan_array[get_local_id(0) + SCAN_BLOCKSIZE];
    if (get_local_id(0) == 0) {
        if (zeroff) output[0] = 0;
        if (aux) aux[get_group_id(0)] = scan_array[2 * SCAN_BLOCKSIZE - 1];
    }
}

__kernel void tally_denselist_lengths(int num_lists, int fdense_list_lengths, int fgridcnt, int fgridoff )
{
    uint list = get_group_id(0) * get_local_size(0) + get_local_id(0);
    if ( list >= num_lists ) return;
    int gridTot = fparam.gridTotal;
    bufI(&fbuf, fdense_list_lengths)[list] = bufI(&fbuf, fgridcnt)[(list+1)*gridTot -1] + bufI(&fbuf, fgridoff)[(list+1)*gridTot -1];
}

__kernel void countingSortFull(int pnum) {
    uint i = get_global_id(0);
    if (i >= pnum) return;
    if (fparam.debug > 1 && i == 0) printf("\ncountingSortFull(): pnum=%u\n", pnum);
    uint icell = bufI(&ftemp, FGCELL)[i];
    if (icell != GRID_UNDEF) {
        uint indx = bufI(&ftemp, FGNDX)[i];
        int sort_ndx = bufI(&fbuf, FGRIDOFF)[icell] + indx;
        float3 zero = (float3)(0, 0, 0);

        bufI(&fbuf, FGRID)[sort_ndx] = sort_ndx;
        bufF3(&fbuf, FPOS)[sort_ndx] = bufF3(&ftemp, FPOS)[i];
        bufF3(&fbuf, FVEL)[sort_ndx] = bufF3(&ftemp, FVEL)[i];
        bufF3(&fbuf, FVEVAL)[sort_ndx] = bufF3(&ftemp, FVEVAL)[i];
        bufF3(&fbuf, FFORCE)[sort_ndx] = zero;
        bufF(&fbuf, FPRESS)[sort_ndx] = bufF(&ftemp, FPRESS)[i];
        bufF(&fbuf, FDENSITY)[sort_ndx] = bufF(&ftemp, FDENSITY)[i];
        bufI(&fbuf, FAGE)[sort_ndx] = bufI(&ftemp, FAGE)[i];
        bufI(&fbuf, FCLR)[sort_ndx] = bufI(&ftemp, FCLR)[i];
        bufI(&fbuf, FGCELL)[sort_ndx] = icell;
        bufI(&fbuf, FGNDX)[sort_ndx] = indx;
        float3 pos = bufF3(&ftemp, FPOS)[i];
        for (int a = 0; a < BONDS_PER_PARTICLE; a++) {
            uint j = bufI(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND];
            uint j_sort_ndx = UINT_MAX;
            uint jcell = GRID_UNDEF;
            if (j < pnum) {
                jcell = bufI(&ftemp, FGCELL)[j];
                uint jndx = UINT_MAX;
                if (jcell != GRID_UNDEF) {
                    jndx = bufI(&ftemp, FGNDX)[j];
                    if ((bufI(&fbuf, FGRIDOFF)[jcell] + jndx) < pnum) {
                        j_sort_ndx = bufI(&fbuf, FGRIDOFF)[jcell] + jndx;
                    }
                }
            }
            bufI(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND] = j_sort_ndx;
            for (int b = 1; b < 5 / DATA_PER_BOND; b++) {
                bufF(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND + b] =
                bufF(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + b];
            }
            bufI(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND + 5] =
            bufI(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 5];
            bufI(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND + 6] =
            bufI(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 6];
            bufF(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND + 7] =
            bufF(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 7];
            bufI(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND + 8] =
            bufF(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 8];
            bufI(&fbuf, FELASTIDX)[sort_ndx * BOND_DATA + a * DATA_PER_BOND + 8] =
            bufI(&ftemp, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 8];

        }
            for (int a = 0; a < BONDS_PER_PARTICLE; a++) {
                uint k = bufI(&ftemp, FPARTICLEIDX)[i * BONDS_PER_PARTICLE * 2 + a * 2];
                uint b = bufI(&ftemp, FPARTICLEIDX)[i * BONDS_PER_PARTICLE * 2 + a * 2 + 1];
                uint kndx, kcell, ksort_ndx = UINT_MAX;
                if (k < pnum) {
                    kcell = bufI(&ftemp, FGCELL)[k];
                    if (kcell != GRID_UNDEF) {
                        kndx = bufI(&ftemp, FGNDX)[k];
                        ksort_ndx = bufI(&fbuf, FGRIDOFF)[kcell] + kndx;
                    }
                }
                bufI(&fbuf, FPARTICLEIDX)[sort_ndx * BONDS_PER_PARTICLE * 2 + a * 2] = ksort_ndx;
                bufI(&fbuf, FPARTICLEIDX)[sort_ndx * BONDS_PER_PARTICLE * 2 + a * 2 + 1] = b;
                bufI(&ftemp, FPARTICLEIDX)[i * BONDS_PER_PARTICLE * 2 + a * 2] = UINT_MAX;
            }

            bufI(&fbuf, FPARTICLE_ID)[sort_ndx] = bufI(&ftemp, FPARTICLE_ID)[i];
            bufI(&fbuf, FMASS_RADIUS)[sort_ndx] = bufI(&ftemp, FMASS_RADIUS)[i];
            bufI(&fbuf, FNERVEIDX)[sort_ndx] = bufI(&ftemp, FNERVEIDX)[i];

            uint* fbuf_epigen = &bufI(&fbuf, FEPIGEN)[sort_ndx];
            uint* ftemp_epigen = &bufI(&ftemp, FEPIGEN)[i];
            for (int a = 0; a < NUM_GENES; a++) fbuf_epigen[pnum * a] = ftemp_epigen[pnum * a];

            float* fbuf_conc = &bufF(&fbuf, FCONC)[sort_ndx * NUM_TF];
            float* ftemp_conc = &bufF(&ftemp, FCONC)[i * NUM_TF];
            for (int a = 0; a < NUM_TF; a++) fbuf_conc[a] = ftemp_conc[a];
        }
}
