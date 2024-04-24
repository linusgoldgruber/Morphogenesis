
| 1. SETUP PIELINE: set the name of your "demo" folder within the Morphogenesis/src/ directory in the config.json file (line 4)


| 2. EXECUTION PIPELINE: createDemo -> checkDemo -> createDemo2

    NB: Creating multiple executables in CMakeLists.txt runs into issues with VTK atm,
        thus example CMakeLists.txt files are provided (see Morphogenesis/src/ direcotry).


| 3. LAUNCH SETTINGS:


createDemo:
-----------------------------------------------
Launch-Configurations:

    Executable:
        build/createDemo

    Arguments:          num_particles, spacing, x_dim, y_dim, z_dim, demoType, simSpace, JSON file

        for example:         20            1       4      4      4       0       6       config.json      (Linus' basic demo)
                 or:        125            1       6      6      6       0       5       config.json      ("free falling" from NH89's github)

----------------------------------------------

checkDemo:
----------------------------------------------
Launch-Configurations:

    Executable:
        build/checkDemo

    Arguments:          JSON file

        for example:   config.json

----------------------------------------------


createDemo2:
----------------------------------------------
Launch-Configurations:

    Executable:
        build/createDemo2

    Arguments:          JSON file

        for example:   config.json

----------------------------------------------


| 4. Environment variables:
----------------------------------------------
COLLECT_GCC_OPTIONS=            '-E' '-v' '-o' '/dev/null' '-mtune=generic' '-march=x86-64'

COMPILER_PATH=                  /usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/:/usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/

CPATH=                          /usr/local/cuda/include:/usr/local/include:/usr/include/x86_64-linux-gnu/c++/11:/usr/include/c++/11:

C_INCLUDE_PATH=                 /usr/local/cuda/include:/usr/local/include:/usr/include:

C_PLUS_INCLUDE_PATH=            /usr/local/cuda/include:/usr/local/include:/usr/include/x86_64-linux-gnu/c++/11:/usr/include/c++/11/tr1:/usr/include/c++/11:

LD_LIBRARY_PATH=$LIBRARY_PATH:  /opt/rocm/opencl/lib:/opt/intel/oneapi/lib

LIBRARY_PATH=                   /usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/11/../../../x86_64-linux-gnu/:/usr/lib/gcc/x86_64-linux-gnu/11/../../../../lib/:/lib/x86_64-linux-gnu/:/lib/../lib/:/usr/lib/x86_64-linux-gnu/:/usr/lib/../lib/:/usr/lib/gcc/x86_64-linux-gnu/11/../../../:/lib/:/usr/lib/

-----------------------------------------------











NICK (old)
0) Allocate buffers (might be automatic, check)
    (i)  particles   FELASTIDX, FNERVEIDX, FCONC, FEPIGEN
    (ii) genome
    

1) Initalize correct UID for each particle

2) Initialize Buffer[FELASTIDX],   FNERVEIDX   , FCONC , FEPIGEN

3) UpdateGenome(); //  sends genome to device. // NB need to initialize genome from file, or something.

4) Need script to generate simulations + functions to read them from file. ? what are the available i/o functions ? currently FluidSystem::SavePoints* () 
Also see how files are read in OpenCL-SPH, and how I planned to use Json, Yaml, hdf5 etc...
Probably use a folder with named csv files.

5) 


-----------------------------------------------------------------------------------------------------------------------------------------------------
CMAKE:

Where is the source code: 
/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src

Where to build the binaries:
/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/build

----
load_sim.cpp

Executable:
/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/build/load_sim

Arguments:
/home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/demo /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/demo/out 3 10 n n y y 3 n n /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/opencl_config.json


| Code Repository--------------------------------------------------------------------------------------------------------------------------


INSERT PARTICLES KERNEL:


__kernel void insertParticlesCL(
    __global struct FParams* m_FParamsDevice,
//     __global struct FBufs* m_FluidDevice,
    int pnum,
    volatile __global float4* fpos,
    volatile __global uint* fgcell,
    volatile __global uint* fgndx,
    volatile __global uint* fepigen,
    volatile __global int* fgridcnt,
    volatile __global int* fgridcnt_active_genes
    )
{
uint i = get_global_id(0);
    if ( i >= pnum ) return;

    float4  gridMin     =  m_FParamsDevice->gridMin;
    float4  gridDelta   =  m_FParamsDevice->gridDelta;
    int4    gridRes     =  m_FParamsDevice->gridRes;
    int4    gridScan    =  m_FParamsDevice->gridScanMax;
    int     gridTot     =  m_FParamsDevice->gridTotal;

    float4 gcf = (fpos[i] - gridMin) * gridDelta;

//     printf("GCF after Calculation: Thread ID: %u, gcf: (%f, %f, %f), gridDelta: (%f, %f, %f)\n",
//         i, gcf.x, gcf.y, gcf.z, gridDelta.x, gridDelta.y, gridDelta.z);
    int4 gc  = (int4)(gcf.x, gcf.y, gcf.z, gcf.w);
    int gs  = (gc.y * gridRes.z + gc.z) * gridRes.x + gc.x;

//     printf("Thread ID: %u, m_FParamsDevice->gridDelta: (%f, %f, %f)\n",
//         i, m_FParamsDevice->gridDelta.x, m_FParamsDevice->gridDelta.y, m_FParamsDevice->gridDelta.z);

//     printf("Thread ID: %u, float4 fpos[%u]: (%v4hlf)\n",
//         i, i, fpos[i]);

// printf("Thread ID: %u, m_FParamsDevice->gridTotal: (%u)\n",
//        i, m_FParamsDevice->gridTotal);

    if(i==pnum-1) printf("\n\ninsertParticles()1: gridTot=%i,  i=%u: gc.x=%i, gc.y=%i, gc.z=%i, gs=%i \t gridScan.x=%i, gridScan.y=%i, gridScan.z=%i, gridTot=%u,\t gridDelta=(%f,%f,%f) gridMin=(%f,%f,%f) gridRes=(%i,%i,%i)",
    gridTot, i, gc.x, gc.y, gc.z, gs,  gridScan.x, gridScan.y, gridScan.z, gridTot, gridDelta.x, gridDelta.y, gridDelta.z,  gridMin.x, gridMin.y, gridMin.z, gridRes.x, gridRes.y, gridRes.z );

    if ( gc.x >= 1 && gc.x <= gridScan.x && gc.y >= 1 && gc.y <= gridScan.y && gc.z >= 1 && gc.z <= gridScan.z ) {

		fgcell[i] = gs;                                    // Grid cell insert.
// 		printf("InsideLoop: Thread ID: %u, fgcell: %u\n",             // all zero (0)
//         i, fgcell[i]);
//         printf("Thread ID: %u, fgridcnt[%d] before atomic_add: %d\n", i, gs, fgridcnt[gs]);
        fgndx[i] = atomic_add(&fgridcnt[gs], 1);          // Grid counts.
//         printf("Thread ID: %u, fgridcnt[%d] after atomic_add: %d\n", i, gs, fgridcnt[gs]);
//         printf("Thread ID: %u, fgndx[%d]: %u\n", i, i, fgndx[i]);
        //  ## add counters for dense lists. ##############
        // for each gene, if active, then atomicAdd bin count for gene
        for(int gene=0; gene<NUM_GENES; gene++){ // NB data ordered FEPIGEN[gene][particle] AND +ve int values -> active genes.
            //if(m_FParamsDevice->debug>2 && i==0)printf("\n");
            if (fepigen[i + gene * m_FParamsDevice->maxPoints] > 0) {  // "if((int)bufI(m_FluidDevice, FEPIGEN)" may clash with INT_MAX
                atomic_add( &fgridcnt_active_genes[gene*gridTot +gs], 1 );
                //if(fparam.debug>2 && (gene==6||gene==9) /*i<10*/) printf("\ninsertParticles()2: i=,%u, gene=,%u, gs=,%u, fbuf.bufI(FGRIDCNT_ACTIVE_GENES)[ gene*gridTot  + gs ]=,%u",
                //    i, gene, gs, fbuf.bufI(FGRIDCNT_ACTIVE_GENES)[ gene*gridTot  + gs ]);
            }
            // could use a small array of uints to store gene activity as bits. This would reduce the reads, but require bitshift and mask to read.
            //if(fparam.debug>2 && i==0)printf("\ninsertParticles()3: fbuf.bufI(FEPIGEN) [i*NUM_GENES + gene]=%u  gene=%u  i=%u,",fbuf.bufI(FEPIGEN)[gene*pnum + i/* i*NUM_GENES + gene*/], gene ,i  );
        }
    } else {
        fgcell[i] = GRID_UNDEF;
//         printf("Thread ID: %u, GRID_UNDEF\n",i);

    }

