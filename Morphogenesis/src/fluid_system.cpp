#include <assert.h>
#include <iostream>
#include <CL/cl.h>
#include <stdlib.h>
#include <unistd.h>
#include <chrono>
#include <cstring>
#include "fluid_system.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
#include "RunCL.h"
#include "fluid.h"
#include "fluid_system.h"
#include "RunCL.cpp"
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }
//#include <curand_kernel.h> //NOT REPLACED YET. POSSIBLE OPTIONS: https://acesse.dev/Y9Zcq

using namespace std;

string checkerror(int input) {
		int errorCode = input;
		switch (errorCode) {
		case -9999:											return "Illegal read or write to a buffer";		// NVidia error code
		case CL_DEVICE_NOT_FOUND:							return "CL_DEVICE_NOT_FOUND";
		case CL_DEVICE_NOT_AVAILABLE:						return "CL_DEVICE_NOT_AVAILABLE";
		case CL_COMPILER_NOT_AVAILABLE:						return "CL_COMPILER_NOT_AVAILABLE";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:				return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case CL_OUT_OF_RESOURCES:							return "CL_OUT_OF_RESOURCES";
		case CL_OUT_OF_HOST_MEMORY:							return "CL_OUT_OF_HOST_MEMORY";
		case CL_PROFILING_INFO_NOT_AVAILABLE:				return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case CL_MEM_COPY_OVERLAP:							return "CL_MEM_COPY_OVERLAP";
		case CL_IMAGE_FORMAT_MISMATCH:						return "CL_IMAGE_FORMAT_MISMATCH";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:					return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case CL_BUILD_PROGRAM_FAILURE:						return "CL_BUILD_PROGRAM_FAILURE";
		case CL_MAP_FAILURE:								return "CL_MAP_FAILURE";
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:				return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:	return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case CL_INVALID_VALUE:								return "CL_INVALID_VALUE";
		case CL_INVALID_DEVICE_TYPE:						return "CL_INVALID_DEVICE_TYPE";
		case CL_INVALID_PLATFORM:							return "CL_INVALID_PLATFORM";
		case CL_INVALID_DEVICE:								return "CL_INVALID_DEVICE";
		case CL_INVALID_CONTEXT:							return "CL_INVALID_CONTEXT";
		case CL_INVALID_QUEUE_PROPERTIES:					return "CL_INVALID_QUEUE_PROPERTIES";
		case CL_INVALID_COMMAND_QUEUE:						return "CL_INVALID_COMMAND_QUEUE";
		case CL_INVALID_HOST_PTR:							return "CL_INVALID_HOST_PTR";
		case CL_INVALID_MEM_OBJECT:							return "CL_INVALID_MEM_OBJECT";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:			return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case CL_INVALID_IMAGE_SIZE:							return "CL_INVALID_IMAGE_SIZE";
		case CL_INVALID_SAMPLER:							return "CL_INVALID_SAMPLER";
		case CL_INVALID_BINARY:								return "CL_INVALID_BINARY";
		case CL_INVALID_BUILD_OPTIONS:						return "CL_INVALID_BUILD_OPTIONS";
		case CL_INVALID_PROGRAM:							return "CL_INVALID_PROGRAM";
		case CL_INVALID_PROGRAM_EXECUTABLE:					return "CL_INVALID_PROGRAM_EXECUTABLE";
		case CL_INVALID_KERNEL_NAME:						return "CL_INVALID_KERNEL_NAME";
		case CL_INVALID_KERNEL_DEFINITION:					return "CL_INVALID_KERNEL_DEFINITION";
		case CL_INVALID_KERNEL:								return "CL_INVALID_KERNEL";
		case CL_INVALID_ARG_INDEX:							return "CL_INVALID_ARG_INDEX";
		case CL_INVALID_ARG_VALUE:							return "CL_INVALID_ARG_VALUE";
		case CL_INVALID_ARG_SIZE:							return "CL_INVALID_ARG_SIZE";
		case CL_INVALID_KERNEL_ARGS:						return "CL_INVALID_KERNEL_ARGS";
		case CL_INVALID_WORK_DIMENSION:						return "CL_INVALID_WORK_DIMENSION";
		case CL_INVALID_WORK_GROUP_SIZE:					return "CL_INVALID_WORK_GROUP_SIZE";
		case CL_INVALID_WORK_ITEM_SIZE:						return "CL_INVALID_WORK_ITEM_SIZE";
		case CL_INVALID_GLOBAL_OFFSET:						return "CL_INVALID_GLOBAL_OFFSET";
		case CL_INVALID_EVENT_WAIT_LIST:					return "CL_INVALID_EVENT_WAIT_LIST";
		case CL_INVALID_EVENT:								return "CL_INVALID_EVENT";
		case CL_INVALID_OPERATION:							return "CL_INVALID_OPERATION";
		case CL_INVALID_GL_OBJECT:							return "CL_INVALID_GL_OBJECT";
		case CL_INVALID_BUFFER_SIZE:						return "CL_INVALID_BUFFER_SIZE";
		case CL_INVALID_MIP_LEVEL:							return "CL_INVALID_MIP_LEVEL";
		case CL_INVALID_GLOBAL_WORK_SIZE:					return "CL_INVALID_GLOBAL_WORK_SIZE";
        #if CL_HPP_MINIMUM_OPENCL_VERSION >= 200
            case CL_INVALID_DEVICE_QUEUE:						return "CL_INVALID_DEVICE_QUEUE";
            case CL_INVALID_PIPE_SIZE:							return "CL_INVALID_PIPE_SIZE";
        #endif
            default:											return "unknown error code";
		}
	}

bool clCheck(cl_int status, const char* method, const char* apicall, const char* arg, bool bDebug)
{

    // DEBUG IMPLEMENTATION MISSING!!!

    if (status != CL_SUCCESS) {
        std::string errorMessage = checkerror(status);
        std::cout << "OpenCL Error: " << errorMessage << std::endl;
        std::cout << "Caller: " << method << std::endl;
        std::cout << "Call: " << apicall << std::endl;
        std::cout << "Args: " << arg << std::endl;
        return false;
    }
    return true;
}
/*
FluidSystem::FluidSystem (RunCL& runcl){
    if (m_FParams.debug>1)cout<<"\n\nFluidSystem ()"<<std::flush;
    memset ( &m_Fluid, 0,		sizeof(FBufs) );
    memset ( &m_FluidTemp, 0,	sizeof(FBufs) );
    memset ( &m_FParams, 0,		sizeof(FParams) );
    memset ( &m_FGenome, 0,		sizeof(FGenome) );
    mNumPoints = 0;
    mMaxPoints = 0;
    mPackGrid = 0x0;
    m_Frame = 0;
    //Kernels are get defined in RunCl.h!!!
    //for (int n=0; n < FUNC_MAX; n++ ) runcl.m_Kern[n] = (cl_kernel) -1;
}
*/
bool FluidSystem::clCheck(cl_int status, const char* method, const char* apicall, const char* arg, bool bDebug)
{
    // DEBUG IMPLEMENTATION MISSING!!!

    if (status != CL_SUCCESS) {
        std::string errorMessage = checkerror(status);
        std::cout << "OpenCL Error: " << errorMessage << std::endl;
        std::cout << "Caller: " << method << std::endl;
        std::cout << "Call: " << apicall << std::endl;
        std::cout << "Args: " << arg << std::endl;
        return false;
    }
    return true;
}

/*void FluidSystem::LoadKernel ( int fid, std::string func ){
    char cfn[512];
    strcpy ( cfn, func.c_str() );

    if ( m_Kern[fid] == (cl_kernel) -1 )
        clCheck ( cuModuleGetFunction ( &m_Kern[fid], m_Program, cfn ), "LoadKernel", "cuModuleGetFunction", cfn, mbDebug );
}*/

//Kernels already get defined in RunCL.h!!!
/*
void FluidSystem::LoadKernel(RunCL& runcl, int fid, std::string func) {
    if (runcl.m_Kern[fid] == nullptr) {
        cl_int err;
        runcl.m_Kern[fid] = clCreateKernel(m_program, func.c_str(), &err);
        clCheck(err == CL_SUCCESS, "LoadKernel", "clCreateKernel", func.c_str(), mbDebug);
    }
}*/


void FluidSystem::Initialize(RunCL& runcl){             //Left aside for now, implement by copying from InitializeOpenCLused for CPU only for "check_demo".

    if (m_FParams.debug>1)std::cout << "FluidSystem::Initialize() \n";
    // An FBufs struct holds an array of pointers.
    // Clear all buffers
    memset ( &m_Fluid, 0,		sizeof(FBufs) );
    memset ( &m_FluidTemp, 0,	sizeof(FBufs) );
    memset ( &m_FParams, 0,		sizeof(FParams) );
    memset ( &m_FGenome, 0,		sizeof(FGenome) );
    mNumPoints = 0;
    mMaxPoints = 0;
    mPackGrid = 0x0;
    m_Frame = 0;
    cout << "-----Initialization successful-----" << endl;

    if (m_FParams.debug>1)std::cout << "Chk1.4 \n";

    // Allocate the sim parameters CUCLCUCL
    //AllocateBuffer ( FPARAMS,		sizeof(FParams),	0,	1,	 GPU_OFF,     CPU_YES );//AllocateBuffer ( int buf_id, int stride,     int cpucnt, int gpucnt,    int gpumode,    int cpumode )

    if (m_FParams.debug>1)std::cout << "Chk1.5 \n";

    m_Time = 0;
    mNumPoints = 0;			// reset count

    if (m_FParams.debug>1)std::cout << "Chk1.6 \n";
}

/* STILL NECESSARY??????
void FluidSystem::AllocateBuffer ( int buf_id, int stride, int cpucnt, int gpucnt, int gpumode, int cpumode ){   // mallocs a buffer - called by FluidSystem::Initialize(), AllocateParticles, and AllocateGrid()
//also called by WriteDemoSimParams(..)
    bool rtn = true;
    if (m_FParams.debug>1)std::cout<<"\nAllocateBuffer ( int buf_id="<<buf_id<<", int stride="<<stride<<", int cpucnt="<<cpucnt<<", int gpucnt="<<gpucnt<<", int "<<gpumode<<", int "<<cpumode<<" )\t"<<std::flush;
    if (cpumode == CPU_YES) {
        char* src_buf  = bufC(&m_Fluid, buf_id);
        cl_mem dest_buf = clCreateBuffer(m_context, CL_MEM_READ_WRITE, cpucnt*stride, NULL, &err); //  ####  malloc the buffer   ####
        if (src_buf != 0x0) {
            clEnqueueWriteBuffer(queue, dest_buf, CL_TRUE, 0, cpucnt*stride, src_buf, 0, NULL, NULL);
            free(src_buf);
        }
        setBuf(&m_Fluid, buf_id, dest_buf); // stores pointer to buffer in mcpu[buf_id]
}*/

void FluidSystem::UpdateGenome (RunCL& runcl){              // Update Genome on GPU
    clCheck ( clEnqueueWriteBuffer(runcl.m_queue, runcl.m_FGenomeDevice, CL_TRUE, 0, sizeof(FGenome), &m_FGenome, 0, NULL, NULL), "FluidUpdateGenome", "clEnqueueWriteBuffer", "m_FGenomeDevice", mbDebug);
}

FGenome	FluidSystem::GetGenome(RunCL& runcl){
            FGenome tempGenome = m_FGenome;
            for (int i=0; i<NUM_GENES; i++)tempGenome.mutability[i]=m_FGenome.mutability[i];
            for (int i=0; i<NUM_GENES; i++)tempGenome.delay[i]=m_FGenome.delay[i];
            for (int i=0; i<NUM_GENES; i++)for (int j=0; j<NUM_GENES; j++)tempGenome.sensitivity[i][j]=m_FGenome.sensitivity[i][j];

            for (int i=0; i<NUM_TF; i++)tempGenome.tf_diffusability[i]=m_FGenome.tf_diffusability[i];
            for (int i=0; i<NUM_TF; i++)tempGenome.tf_breakdown_rate[i]=m_FGenome.tf_breakdown_rate[i];

            for (int i=0; i<NUM_GENES; i++)for (int j=0; j<2*NUM_TF+1; j++)tempGenome.secrete[i][j]=m_FGenome.secrete[i][j];
            for (int i=0; i<NUM_GENES; i++)for (int j=0; j<2*NUM_GENES+1; j++)tempGenome.activate[i][j]=m_FGenome.activate[i][j];

            for (int i=0; i<3;i++)for(int j=0; j<12; j++)tempGenome.param[i][j]=m_FGenome.param[i][j];
            std::cout<<"\nGetGenome(): m_FGenome.delay[0]="<<m_FGenome.delay[0]<<"\ttempGenome.delay[0]="<<tempGenome.delay[0]<<std::flush;
            return tempGenome;
}

void FluidSystem::FluidParamCL (RunCL& runcl, float ss, float sr, float pr, float mass, float rest, float3 bmin, float3 bmax, float estiff, float istiff, float visc, float surface_tension, float damp, float fmin, float fmax, float ffreq, float gslope, float gx, float gy, float gz, float al, float vl, float a_f, float a_p ){
    m_FParams.psimscale = ss;
    m_FParams.psmoothradius = sr;
    m_FParams.pradius = pr;
    m_FParams.r2 = sr * sr;
    m_FParams.pmass = mass;
    m_FParams.prest_dens = rest;
    m_FParams.pboundmin = bmin;
    m_FParams.pboundmax = bmax;
    m_FParams.pextstiff = estiff;
    m_FParams.pintstiff = istiff;
    m_FParams.pvisc = visc;
    m_FParams.psurface_t = surface_tension;
    m_FParams.pdamp = damp;
    m_FParams.pforce_min = fmin;
    m_FParams.pforce_max = fmax;
    m_FParams.pforce_freq = ffreq;
    m_FParams.pground_slope = gslope;
    m_FParams.pgravity = make_float3( gx, gy, gz );
    m_FParams.AL = al;
    m_FParams.AL2 = al * al;
    m_FParams.VL = vl;
    m_FParams.VL2 = vl * vl;
    //m_FParams.pemit = emit;

    m_FParams.pdist = pow ( m_FParams.pmass / m_FParams.prest_dens, 1/3.0f );
                                                                                // Normalization constants.
    m_FParams.poly6kern = 315.0f / (64.0f * 3.141592f * pow( sr, 9.0f) );
    m_FParams.wendlandC2kern = 21 / (2 * 3.141592f );   // This is the value calculated in SymPy as per Wendland C2 as per (Dehnen & Aly 2012)
    // 16   // The  WC2 kernel in DualSPHysics assumes  values of 0<=q<=2 , hence the divisor 16pi in the normalisation constant for 3D.
    /* My notes from Sympy my notebook.
    Where Wendland C2 kernel:

        wc2 = (1-r*ss/2*sr)**4  * ((2*q) +1)

    Normalisation constant = 1/integrate( (wc2*(4*pi*r**2)), (r,0, 2*sr/ss)),  NB *(4*pi*r**2) area of a sphere, & 2=basis of wc2.

        =  1/ (288pi - 15552.0πss^2/sr^2 + 77760.0πss^3/sr^3 - 149965.714285714πss^4/sr^4 + 104976.0πss^5/sr^5  )

    */
    /* Notes from DualSPHysics Wiki
    // Normalization const = reciprocal of radial integral of (kernel * area of sphere), found using Sympy.
    // NB using W(r,h)=alpha_D (1-q/2)**4 *(2*q +1), 0<=q<=2, as per DualSPHysics Wiki. Where alpha_D is the normaliation constant.
    // * m_FParams.pmass * m_FParams.psimscale
    */
    m_FParams.spikykern = -45.0f / (3.141592f * pow( sr, 6.0f) );            // spikykern used for force due to pressure.
    m_FParams.lapkern = 45.0f / (3.141592f * pow( sr, 6.0f) );
    // NB Viscosity uses a different kernel, this is the constant portion of its Laplacian.
    // NB Laplacian is a scalar 2nd order differential, "The divergence of the gradient"
    // This Laplacian comes from Muller et al 2003, NB The kernel is defined by the properties of  its Laplacian, gradient and value at the basis (outer limit) of the kernel. The Laplacian is the form used in the code. The equation of the kernel in Muller et al seems to be wrong, but this does not matter.

/*
    // -32*(1 - r)**3 + 12*(1 - r)**2*(4*r + 1)  // the Laplacian of  WC2 = (1-r)**4 *(1+4*r)
//(15*r**2*(h/r**3 + 2/h**2 - 3*r/h**3)/(2*pi*h**3) + 15*r*(-h/(2*r**2) + 2*r/h**2 - 3*r**2/(2*h**3))/(pi*h**3))/r**2
//(45/pi*h^6)((h^2/12r^3)+(2h/3)-(3r/4))

//(r**2*(h/r**3 + 2/h**2 - 3*r/h**3) + 2*r*(-h/(2*r**2) + 2*r/h**2 - 3*r**2/(2*h**3) ) )/r**2
*/

    m_FParams.gausskern = 1.0f / pow(3.141592f * 2.0f*sr*sr, 3.0f/2.0f);     // Gaussian not currently used.

    m_FParams.H = m_FParams.psmoothradius / m_FParams.psimscale;
    m_FParams.d2 = m_FParams.psimscale * m_FParams.psimscale;
    m_FParams.rd2 = m_FParams.r2 / m_FParams.d2;
    m_FParams.vterm = m_FParams.lapkern * m_FParams.pvisc;

    m_FParams.actuation_factor = a_f;
    m_FParams.actuation_period = a_p;


    // Transfer sim params to device
    clCheck ( clEnqueueWriteBuffer(runcl.m_queue, runcl.m_FParamDevice, CL_TRUE, 0, sizeof(FParams), &m_FParams, 0, NULL, NULL), "FluidParamCL", "clEnqueueWriteBuffer", "clFParams", mbDebug);
}


void FluidSystem::UpdateParams (RunCL& runcl){
    // Update Params on GPU
    Vector3DF grav = Vector3DF_multiplyFloat(&m_Vec[PPLANE_GRAV_DIR], m_Param[PGRAV]);
    /*FluidParamCL (runcl,  m_Param[PSIMSCALE], m_Param[PSMOOTHRADIUS], m_Param[PRADIUS], m_Param[PMASS], m_Param[PRESTDENSITY],
                      *(float3*)& m_Vec[PBOUNDMIN], *(float3*)& m_Vec[PBOUNDMAX], m_Param[PEXTSTIFF], m_Param[PINTSTIFF],
                      m_Param[PVISC], m_Param[PSURFACE_TENSION], m_Param[PEXTDAMP], m_Param[PFORCE_MIN], m_Param[PFORCE_MAX], m_Param[PFORCE_FREQ],
                      m_Param[PGROUND_SLOPE], grav.x, grav.y, grav.z, m_Param[PACCEL_LIMIT], m_Param[PVEL_LIMIT],
                      m_Param[PACTUATION_FACTOR], m_Param[PACTUATION_PERIOD]);*/
}

void FluidSystem::SetParam (RunCL& runcl, int p, float v ){
    m_Param[p] = v;
    UpdateParams (runcl);
}

void FluidSystem::SetVec (RunCL& runcl,  int p, Vector3DF v ){
    m_Vec[p] = v;
    UpdateParams (runcl);
}
/*      Probably unnecessary, no desctroctor needed in OpenCL, rest is done in RunCL::CleanUp()
void FluidSystem::Exit (RunCL& runcl){
    // Free fluid buffers
    clCheck(clFinish(runcl.m_queue), "Exit ", "clFinish", "before cudaDeviceReset()", mbDebug);
    for (int n=0; n < MAX_BUF; n++ ) {
        if (m_FParams.debug>0)std::cout << "\n n = " << n << std::flush;
        if ( bufC(&m_fluid, n) != 0x0 )
            free ( bufC(&m_fluid, n) );
    }
    //size_t   free1, free2, total;
    //cudaMemGetInfo(&free1, &total);
    cl_ulong free1, free2, total;
    clGetDeviceInfo(devices[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(total), &total, NULL);

    // Free memory cant be checked so easily in OpenCl, need to implement the funtionailty manually, e.g. by a varialbe called AllocatedMem CUCLCUCL
    if (m_FParams.debug>0)printf("\nCuda Memory, before cudaDeviceReset(): free=%lu, total=%lu.\t",free1,total);
    //clCheck(clFinish(), "Exit ", "clFinish", "before cudaDeviceReset()", mbDebug);
    if(program != 0x0){
        if (m_FParams.debug>0)printf("\nclRelease()\n");
        //cudaDeviceReset(); // Destroy all allocations and reset all state on the current device in the current process. // must only operate if we have a cuda instance.
        clReleaseMemObject(clFBuf);
        clReleaseMemObject(clFTemp);
        clReleaseMemObject(clFParams);
        clReleaseMemObject(clFGenome);
        clReleaseProgram(prograudaDeviceResetm);
        clReleaseKernel(kernel);
        clReleaseCommandQueue(queue);
        clReleaseContext(clContext);
    }

    cudaMemGetInfo(&free2, &total);
    if (m_FParams.debug>0)printf("\nAfter cudaDeviceReset(): free=%lu, total=%lu, released=%lu.\n",free2,total,(free2-free1) );
    exit(0);
}*/
