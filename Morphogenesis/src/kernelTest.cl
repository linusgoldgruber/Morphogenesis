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
        .freezeBoolToInt = 0,
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

__constant struct FPrefix       fprefix = {};
//__constant_FBondParams    fbondparams;    // GPU copy of remodelling parameters.
//__constant uint			gridActive;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// Kernels /////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__kernel void memset32d_kernel(__global int* buffer, const int value) {
    int gid = get_global_id(0);
    buffer[gid] = value;
}

__kernel void prefixFixup(__global uint *input, __global uint *aux, int len) {

    unsigned int t = get_local_id(0);
    unsigned int start = t + 2 * get_group_id(0) * SCAN_BLOCKSIZE;
    if (start < len)					input[start] += aux[get_group_id(0)];
    if (start + SCAN_BLOCKSIZE < len)   input[start + SCAN_BLOCKSIZE] += aux[get_group_id(0)];
}

__kernel void prefixSum(
        __global uint *input,
        __global uint *output,
        __global uint *aux,
        __global size_t *offset_array0,
        __global size_t *offset_scan0,
        int len,
        int zeroff
        ) //, __constant size_t *offset_scan0)
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

float contributePressure (int i, float3 p, int cell, float *sum_p6k) {

    if ( bufI(&fbuf, FGRIDCNT)[cell] == 0 ) return 0.0;                       // If the cell is empty, skip it.

    float3 dist;
    float dsq, r, q, b, c, sum = 0.0;
    float d2 = fparam.psimscale * fparam.psimscale;
    float r2 = fparam.r2;
    float sr = fparam.psmoothradius;

    int clast = bufI(&fbuf, FGRIDOFF)[cell] + bufI(&fbuf, FGRIDCNT)[cell];

    for (int cndx = bufI(&fbuf, FGRIDOFF)[cell]; cndx < clast; cndx++) {
        int pndx = bufI(&fbuf, FGRID)[cndx];
        dist = p - bufF3(&fbuf, FPOS)[pndx];
        dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);

        if (dsq < r2 && dsq > 0.0) {
            r = sqrt(dsq);
            q = r / sr;
            b = (1 - q);
            b *= b;
            b *= b;
            sum += b * (4 * q + 1);

            c = (r2 - dsq) * d2;
            sum_p6k[i] += c * c * c;
        }
    }
    return sum;;
}

__kernel void computePressure ( int pnum ) {
    int i = get_global_id(0); // particle index
    if (i >= pnum) return;

    // Get search cell
    int nadj = (1 * fparam.gridRes.z + 1) * fparam.gridRes.x + 1;
    uint gc = bufI(&fbuf, FGCELL)[i]; // get grid cell of the current particle.
    if (gc == GRID_UNDEF) return;   // IF particle not in the simulation
    gc -= nadj;

    // Sum Pressures
    float3 pos = bufF3(&fbuf, FPOS)[i];
    float sum = 0.0;
    float sum_p6k = 0.0;
    for (int c = 0; c < fparam.gridAdjCnt; c++) {
        sum += contributePressure(i, pos, gc + fparam.gridAdj[c], &sum_p6k);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Compute Density & Pressure
    sum = sum * fparam.pmass * fparam.wendlandC2kern;

    if (sum == 0.0) sum = 1.0;
    bufF(&fbuf, FPRESS)[i] = (sum - fparam.prest_dens) * fparam.pintstiff; // pressure = (diff from rest density) * stiffness
    bufF(&fbuf, FDENSITY)[i] = 1.0f / sum;
}

float3 contributeForce ( int i, float3 ipos, float3 iveleval, float ipress, float idens, int cell)
{
	if ( bufI(&fbuf, FGRIDCNT)[cell] == 0 ) return (float3)(0,0,0);                                        // If the cell is empty, skip it.
	float  dsq, sdist, c, r, sr=fparam.psmoothradius;//1.0;//
    float3 pterm= (float3)(0,0,0), sterm= (float3)(0,0,0), vterm= (float3)(0,0,0), forcej= (float3)(0,0,0), delta_v= (float3)(0,0,0);                                                              // pressure, surface tension and viscosity terms.
	float3 dist     = (float3)(0,0,0),      eterm = (float3)(0,0,0),    force = (float3)(0,0,0);
	uint   j;
	int    clast    = bufI(&fbuf, FGRIDOFF)[cell] + bufI(&fbuf, FGRIDCNT)[cell];                                // index of last particle in this cell
    uint k =0 ;
    for (int cndx = bufI(&fbuf, FGRIDOFF)[cell]; cndx < clast; cndx++ ) {                                     // For particles in this cell.
        k++;
		j           = bufI(&fbuf, FGRID)[ cndx ];
		dist        = ( ipos - bufF3(&fbuf, FPOS)[ j ] );                                                     // dist in cm (Rama's comment)
		dsq         = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);                                      // scalar distance squared
		r           = sqrt(dsq);


        if ( dsq < 1 /*fparam.rd2*/ && dsq > 0) {                                                           // IF in-range && not the same particle
            float kern = pow((sr - r),3);                                                                   // used as a component of surface tension kernel AND directly in viscosity
            sdist   = sqrt(dsq * fparam.d2);                                                                // smoothing distance
            float press = 100*(ipress+bufF(&fbuf, FPRESS)[j]);///sdist

            pterm = idens * bufF(&fbuf, FDENSITY)[j] *  100.0f* (dist/r) *(press*kern - (fparam.psurface_t/*0.4*/)*pow((sr - r),2));       // 1000 = hydroststic stiffness
            delta_v = bufF3(&fbuf, FVEVAL)[j] - iveleval;
            vterm =  100000.0f* delta_v * kern;// (1/2)*pow((sr - r),3) ; // 10000.0 gives fluid, 100000.0 gives visco-elastic behaviour.
            //if (i==1) printf("\n contributeForce : fparam.psurface_t=%f,  fparam.sterm=%f, fparam.pvisc=%f, fparam.vterm=%f ", fparam.psurface_t, fparam.sterm, fparam.pvisc, fparam.vterm );
            /*
             sdist   = sqrt(dsq * fparam.d2);                                                                // smoothing distance = sqrt(dist^2 * sim_scale^2))
             c       = ( fparam.psmoothradius - sdist );
             pterm   = (dist/sdist) * pow((fparam.psmoothradius - sqrt(dsq)), 3) * (fparam.psmoothradius - dsq) ;
             * fparam.psimscale * -0.5f * c * fparam.spikykern   * ( ipress + bufF(&fbuf, FPRESS)[ j ] )/ sdist )  ;       // pressure term
            //sterm   = (dist/dsq) * fparam.sterm * cos(3*CUDART_PI_F*r/(2*fparam.psmoothradius));  // can we use sdist in placeof r ?  or in place od dsq? What about pressure?
			//vterm   =  fparam.vterm * ( bufF3(&fbuf, FVEVAL)[ j ] - iveleval );  // make_float3(0,0,0);//
			forcej  += ( pterm + sterm + vterm) * c * idens * (bufF(&fbuf, FDENSITY)[ j ] );  // fluid force
            */
            force   +=  pterm + vterm  ;
            /*
            if(fparam.debug>0 && i<5 && k<2)  printf("\ncontribForce : debug=%u. i=%u, r=,%f, sr=,%f, (sr-r)^3=,%f, delta_v=,(%f,%f,%f), vterm=(%f,%f,%f), pterm(%f,%f,%f), \t\t press=,%f, sdist=,%f, dsq=,%f, fparam.d2=,%f  kern=,%f, \t\t idens=%f,, bufF(&fbuf, FDENSITY)[j]=,%f, ",fparam.debug, i, r, sr, kern, delta_v.x,delta_v.y,delta_v.z, vterm.x,vterm.y,vterm.z, pterm.x,pterm.y,pterm.z, press, sdist, dsq, fparam.d2, kern, idens, bufF(&fbuf, FDENSITY)[j]);
            */
            /*
            if(i<10) printf("\ncontribForce() : i=,%u, ,cell=,%u,  ,cndx=,%u, ,r=,%f, ,sqrt(fparam.rd2)=r_basis=,%f, ,fparam.psmoothradius=,%f,,sdist=,%f, ,(fparam.psmoothradius-sdist)= c =,%f, \t,ipress=,%f, ,jpress=,%f, ,idens=,%f, ,jdens=,%f,  press=,%f,     \t ,pterm=(,%f,%f,%f,),  ,sterm=(,%f,%f,%f,), ,vterm=(,%f,%f,%f,), ,forcej=(,%f,%f,%f,) ,  ,fparam.vterm=,%f, ,bufF3(&fbuf, FVEVAL)[ j ]=(,%f,%f,%f,), ,iveleval=(,%f,%f,%f,) ",
                i, cell, cndx, r, sqrt(fparam.rd2), fparam.psmoothradius, sdist, c,  ipress, bufF(&fbuf, FPRESS)[j], idens, bufF(&fbuf, FDENSITY)[j], press,   pterm.x,pterm.y,pterm.z, sterm.x,sterm.y,sterm.z, vterm.x,vterm.y,vterm.z, forcej.x,forcej.y,forcej.z,
                fparam.vterm, bufF3(&fbuf, FVEVAL)[j].x, bufF3(&fbuf, FVEVAL)[j].y, bufF3(&fbuf, FVEVAL)[j].z, iveleval.x, iveleval.y, iveleval.z
            );
            */
        }                                                                                                   // end of: IF in-range && not the same particle
    }                                                                                                       // end of loop round particles in this cell
    //if(i<10)  printf("\ncontribForce : i=%u, force=(%f,%f,%f)  ",i, force.x,force.y,force.z  );
    return force;                                                                                           // return fluid force && list of potential bonds fron this cell
}

__kernel void computeForce (int pnum, int freezeBoolToInt, uint frame) {
    uint i = get_global_id(0);
    if (i >= pnum) return;
    uint gc = bufI(&fbuf, FGCELL)[i];
    if (gc == GRID_UNDEF) return;

    gc -= (1 * fparam.gridRes.z + 1) * fparam.gridRes.x + 1;
    float3 force = (float3)(0, 0, 0);
    float3 eterm = (float3)(0, 0, 0);
    float3 dist = (float3)(0, 0, 0);
    float dsq, abs_dist;
    uint bondsToFill = 0;
    uint bonds[BONDS_PER_PARTICLE][2];
    float bond_dsq[BONDS_PER_PARTICLE];
    for (int a = 0; a < BONDS_PER_PARTICLE; a++) {
        bonds[a][0] = UINT_MAX;
        bonds[a][1] = UINT_MAX;
        bond_dsq[a] = fparam.rd2;
    }
    uint i_ID = bufI(&fbuf, FPARTICLE_ID)[i];

    float3 pvel = bufF3(&fbuf, FVEVAL)[i];
    bool hide;
    bool long_bonds = false;
//     for (int a = 0; a < BONDS_PER_PARTICLE; a++) {
//         uint bond = i * BOND_DATA + a * DATA_PER_BOND;
//         uint j = bufI(&fbuf, FELASTIDX)[bond];
//         float restlength = bufF(&fbuf, FELASTIDX)[bond + 2];
//         if (j >= pnum || restlength < 0.000000001) {
//             hide = true;
//             continue;
//         } else hide = false;
//
//         float elastic_limit = bufF(&fbuf, FELASTIDX)[bond + 1];
//         float modulus = bufF(&fbuf, FELASTIDX)[bond + 3];
//         float damping_coeff = bufF(&fbuf, FELASTIDX)[bond + 4];
//         uint other_particle_ID = bufI(&fbuf, FELASTIDX)[bond + 5];
//         uint bondIndex = bufI(&fbuf, FELASTIDX)[bond + 6];
//
//         float3 j_pos = bufF3(&fbuf, FPOS)[j];
//
//         dist = (bufF3(&fbuf, FPOS)[i] - j_pos);
//         dsq = (dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);
//         abs_dist = sqrt(dsq) + FLT_MIN;
//         float3 rel_vel = bufF3(&fbuf, FVEVAL)[j] - pvel;
//
//         float spring_strain = fmax(0.0f, (abs_dist - restlength) / restlength);
//         bufF(&fbuf, FELASTIDX)[bond + /*6*/strain_sq_integrator] = (bufF(&fbuf, FELASTIDX)[bond + /*6*/strain_sq_integrator] + spring_strain * spring_strain);
//         bufF(&fbuf, FELASTIDX)[bond + /*7*/strain_integrator] = (bufF(&fbuf, FELASTIDX)[bond + /*7*/strain_integrator] + spring_strain);
//
//         eterm = ((float)(abs_dist < elastic_limit)) * (((dist / abs_dist) * spring_strain * modulus) - damping_coeff * rel_vel) / (fparam.pmass);
//
//         if (fparam.debug > 0 && abs_dist > 1.5) {
//             long_bonds = true;
//             printf("\ncomputeForce() 1: frame=%u, i=%u, i_ID=%u, j=%u, j_ID=%u, other_particle_ID=%u, bond=%u, eterm=(%f,%f,%f) restlength=%f, modulus=%f , abs_dist=%f , spring_strain=%f , strain_integrator=%f, damping_coeff*rel_vel.z/fparam.pmass=%f, ((dist/abs_dist) * spring_strain * modulus) / fparam.pmass=%f ",
//                 frame, i, i_ID, j, bufI(&fbuf, FPARTICLE_ID)[j], other_particle_ID, a, eterm.x, eterm.y, eterm.z, restlength, modulus, abs_dist, spring_strain, bufF(&fbuf, FELASTIDX)[bond + 7],
//                 damping_coeff * rel_vel.z / fparam.pmass, (((dist.z / abs_dist) * spring_strain * modulus) / fparam.pmass)
//                 );
//         }
//
//         if (isnan(eterm.x) || isnan(eterm.y) || isnan(eterm.z)) {
//             if (!hide) {
//                 printf("\n#### i=%i, j=%i, bond=%i, eterm.x=%f, eterm.y=%f, eterm.z=%f \t####", i,j,a, eterm.x,eterm.y,eterm.z);
//                 printf("\ncomputeForce() chk3: ParticleID=%u, bond=%u, restlength=%f, modulus=%f , abs_dist=%f , spring_strain=%f , strain_integrator=%f  ", bufI(&fbuf, FPARTICLE_ID)[i], a, restlength, modulus, abs_dist, spring_strain, bufF(&fbuf, FELASTIDX)[bond + 7]);
//             }
//         }else {
//             force -= eterm;
//             atomic_fetch_add(bufF3(&fbuf, FFORCE)[j].x += 1.0f);
//             atomic_fetch_add(bufF3(&fbuf, FFORCE)[j].y += 1.0f);
//             atomic_fetch_add(bufF3(&fbuf, FFORCE)[j].z += 1.0f);
//         }
//         if (abs_dist >= elastic_limit && freezeBoolToInt == 0) {
//             bufF(&fbuf, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 2] = 0.0;
//
//             uint bondIndex_ = bufI(&fbuf, FELASTIDX)[i * BOND_DATA + a * DATA_PER_BOND + 6];
//             if (fparam.debug > 2) printf("\n#### Set to broken, i=%i, j=%i, b=%i, bufI(&fbuf, FPARTICLEIDX)[j*BONDS_PER_PARTICLE*2 + b]=UINT_MAX\t####", i,j,bondIndex_);
//             bondsToFill++;
//         }
//
//         barrier(CLK_LOCAL_MEM_FENCE);
//
//     }


    bondsToFill = BONDS_PER_PARTICLE;
    float3 fluid_force_sum = (float3)(0, 0, 0);
    for (int c = 0; c < fparam.gridAdjCnt; c++) {
        float3 fluid_force = (float3)(0, 0, 0);
        fluid_force = contributeForce(i, bufF3(&fbuf, FPOS)[i], bufF3(&fbuf, FVEVAL)[i], bufF(&fbuf, FPRESS)[i], bufF(&fbuf, FDENSITY)[i], gc + fparam.gridAdj[c]);
        fluid_force_sum += fluid_force;
    }
    force += fluid_force_sum * fparam.pmass;
    if (fparam.debug > 0 && long_bonds == true) {
        fluid_force_sum *= fparam.pmass;
        printf("\nComputeForce 2: i=%u, fluid_force_sum=(%f,%f,%f) force=(%f,%f,%f)",
            i, fluid_force_sum.x, fluid_force_sum.y, fluid_force_sum.z, force.x, force.y, force.z);
    }

    bufF3(&fbuf, FFORCE)[i].x += force.x;                                 // atomicAdd req due to other particles contributing forces via incomming bonds.
    bufF3(&fbuf, FFORCE)[i].y += force.y;                                 // NB need to reset FFORCE to zero in  CountingSortFull(..)
    bufF3(&fbuf, FFORCE)[i].z += force.z;                                 // temporary hack, ? better to write a float3 atomicAdd using atomicCAS ?  ########

}

