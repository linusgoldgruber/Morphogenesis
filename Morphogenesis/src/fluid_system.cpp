#include <assert.h>
#include <iostream>
#include <CL/cl.h>
#include <stdlib.h>
#include <unistd.h>
#include <chrono>
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "vector.h"
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
#include "fluid.h"
#include "fluid_system.h"
#include <CL/cl_platform.h>
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }

#define SDK_SUCCESS 0
#define SDK_FAILURE 1
using namespace std;


std::string checkerror(int input) {
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
//         if (bDebug) {
//         kern_stat = clFinish(m_queue);
//     }

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

///////////////////////////////////////////////////////////////// Initialize OpenCL //////////////////////////////////////////////////////////////////////////////////

FluidSystem::FluidSystem(Json::Value obj_)
{
	obj = obj_;
	verbosity = obj["verbosity"].asInt();
	verbosity = 2;
	std::cout << "RunCL::RunCL verbosity = " << verbosity << std::flush;
																						if(verbosity>0) cout << "\nRunCL_chk 0\n" << flush;
	//createFolders( );																	/*Step1: Getting platforms and choose an available one.*/////////
	cl_uint 		numPlatforms;														//the NO. of platforms
	cl_platform_id 	platform 		= NULL;												//the chosen platform
	cl_int			status 			= clGetPlatformIDs(0, NULL, &numPlatforms);			if (status != CL_SUCCESS){ cout << "Error: Getting platforms!" << endl; exit_(status); }
	uint			conf_platform	= obj["opencl_platform"].asUInt();					if(verbosity>0) cout << "numPlatforms = " << numPlatforms << "\n" << flush;
	if (numPlatforms > conf_platform){																/*Choose the platform.*/
		cl_platform_id* platforms 	= (cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
		status 	 					= clGetPlatformIDs(numPlatforms, platforms, NULL);	if (status != CL_SUCCESS){ cout << "Error: Getting platformsIDs" << endl; exit_(status); }
		platform 					= platforms[ conf_platform ];
		free(platforms);																if(verbosity>0) cout << "\nplatforms[0] = "<<platforms[0]<<", \nplatforms[1] = "<<platforms[1]\
																						<<"\nSelected platform number :"<<conf_platform<<", cl_platform_id platform = " << platform<<"\n"<<flush;
	} else {cout<<"Platform num "<<conf_platform<<" not available."<<flush; exit(0);}

	cl_uint				numDevices = 0;													/*Step 2:Query the platform.*//////////////////////////////////
	cl_device_id        *devices;
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);		if (status != CL_SUCCESS) {cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	uint conf_device = obj["opencl_device"].asUInt();

	if (numDevices > conf_device){
		devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
		status  = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
	}																					if (status != CL_SUCCESS) {cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}

																						if(verbosity>0) cout << "RunCL_chk 3\n" << flush; cout << "cl_device_id  devices = " << devices << "\n" << flush;
	cl_context_properties cps[3]={CL_CONTEXT_PLATFORM,(cl_context_properties)platform,0};/*Step 3: Create context.*////////////////////////////////////
	m_context = clCreateContextFromType( cps, CL_DEVICE_TYPE_GPU, NULL, NULL, &status); if(status!=0) 			{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);}
	deviceId  = devices[conf_device];													/*Step 4: Create command queue & associate context.*///////////
	cl_command_queue_properties prop[] = { 0 };											//  NB Device (GPU) queues are out-of-order execution -> need synchronization.
	m_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);		if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);};
	uload_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);	if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);};
	dload_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);	if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);};
	track_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);	if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);};
																						if(verbosity>0) cout << "RunCL_chk 4: \n" << uload_queue << flush;

																						// Multiple queues for latency hiding: Upload, Download, Mapping, Tracking,... autocalibration, SIRFS, SPMP
																						// NB Might want to create command queues on multiple platforms & devices.
																						// NB might want to divde a task across multiple MPI Ranks on a multi-GPU WS or cluster.
	const char *filename = obj["kernel_filepath"].asCString();							/*Step 5: Create program object*///////////////////////////////
	string sourceStr;
	status 						= convertToString(filename, sourceStr);					if(status!=CL_SUCCESS)	{cout<<"\nconvertToString status="<<checkerror(status)<<"\n"<<flush;exit_(status);}
	const char 	*source 		= sourceStr.c_str();
	size_t 		sourceSize[] 	= { strlen(source) };
	m_program 	= clCreateProgramWithSource(m_context, 1, &source, sourceSize, NULL);

	const char includeOptions[] = "-I /lib/i386-linux-gnu -I \"/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src\""; // Paths to include directories ///// -I /usr/lib/gcc/x86_64-linux-gnu/11/include /////

	status = clBuildProgram(m_program, 1, devices, includeOptions, NULL, NULL);					/*Step 6: Build program.*/////////////////////
	if (status != CL_SUCCESS){
		printf("\nclBuildProgram failed: %d\n", status);
		char buf[0x10000];
		clGetProgramBuildInfo(m_program, deviceId, CL_PROGRAM_BUILD_LOG, 0x10000, buf, NULL);
		printf("\n%s\n End of clBuildProgram error log..", buf);
		exit_(status);
	}

	//Initialize the Array of Kernels "m_Kern[]"
	InitializeKernels(m_program);

	m_FParamDevice=m_FluidDevice=m_FluidTempDevice=m_FGenomeDevice=0;		// set device pointers to zero
																						if(verbosity>0) cout << "\nRunCL_constructor finished\n" << flush;
}

FluidSystem::~FluidSystem()
{
	cl_int status;

//	status = clReleaseKernel(m_Kern[FUNC_INSERT]); 							if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_INSERT status = " << checkerror(status) << "\n" << flush;}
	status = clReleaseKernel(m_Kern[FUNC_COUNTING_SORT]); 					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COUNTING_SORT status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_QUERY]); 							if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_QUERY status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_PRESS]); 					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_PRESS status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_FORCE]); 					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_FORCE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_ADVANCE]); 						if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_ADVANCE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_EMIT]); 							if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_EMIT status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_RANDOMIZE]); 						if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_RANDOMIZE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_SAMPLE]); 							if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_SAMPLE status = " << checkerror(status) << "\n" << flush;}
	status = clReleaseKernel(m_Kern[FUNC_FPREFIXSUM]); 						if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_FPREFIXSUM status = " << checkerror(status) << "\n" << flush;}
	status = clReleaseKernel(m_Kern[FUNC_FPREFIXUP]); 					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_FPREFIXFIXUP status = " << checkerror(status) << "\n" << flush;}
	status = clReleaseKernel(m_Kern[FUNC_TALLYLISTS]);						if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_TALLYLISTS status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_DIFFUSION]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_DIFFUSION status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COUNT_SORT_LISTS]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COUNT_SORT_LISTS status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_GENE_ACTION]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_GENE_ACTION status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_TALLY_GENE_ACTION]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_TALLY_GENE_ACTION status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_BOND_CHANGES]);			if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_BOND_CHANGES status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COUNTING_SORT_CHANGES]);			if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COUNTING_SORT_CHANGES status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_NERVE_ACTION]);			if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_NERVE_ACTION status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_COMPUTE_MUSCLE_CONTRACTION]);		if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_COMPUTE_MUSCLE_CONTRACTION status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_CLEAN_BONDS]);						if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_CLEAN_BONDS status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_HEAL]);							if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_HEAL status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_LENGTHEN_MUSCLE]);					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_LENGTHEN_MUSCLE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_LENGTHEN_TISSUE]);					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_LENGTHEN_TISSUE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_SHORTEN_MUSCLE]);					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_SHORTEN_MUSCLE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_SHORTEN_TISSUE]);					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_SHORTEN_TISSUE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_STRENGTHEN_TISSUE]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_STRENGTHEN_TISSUE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_WEAKEN_MUSCLE]);					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_WEAKEN_MUSCLE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_WEAKEN_TISSUE]);					if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_WEAKEN_TISSUE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_EXTERNAL_ACTUATION]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_EXTERNAL_ACTUATION status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_FIXED]);							if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_FIXED status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_INIT_FCURAND_STATE]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_INIT_FCURAND_STATE status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_ASSEMBLE_MUSCLE_FIBRES_OUTGOING]);	if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_ASSEMBLE_MUSCLE_FIBRES_OUTGOING status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_ASSEMBLE_MUSCLE_FIBRES_INCOMING]);	if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_ASSEMBLE_MUSCLE_FIBRES_INCOMING status = " << checkerror(status) << "\n" << flush;}
// 	status = clReleaseKernel(m_Kern[FUNC_INITIALIZE_BONDS]);				if (status != CL_SUCCESS) {cout << "\nRelease Kernel FUNC_INITIALIZE_BONDS status = " << checkerror(status) << "\n" << flush;}

	status = clReleaseProgram(m_program);			if (status != CL_SUCCESS)	{ cout << "\nRelease Program status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(m_queue);		if (status != CL_SUCCESS)	{ cout << "\nRelease CQ1 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(uload_queue);	if (status != CL_SUCCESS)	{ cout << "\nRelease CQ2 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(dload_queue);	if (status != CL_SUCCESS)	{ cout << "\nRelease CQ3 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(track_queue);	if (status != CL_SUCCESS)	{ cout << "\nRelease CQ4 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseContext(m_context);			if (status != CL_SUCCESS)	{ cout << "\nRelease Context status = " << checkerror(status) <<"\n"<<flush; }
}

void FluidSystem::InitializeOpenCL()
{
																				if(verbosity>0) cout << "RunCL::allocatemem_chk0" << flush;
	stringstream 		ss;
	ss << "allocatemem";
	cl_int status;
	cl_event writeEvt;

	int layerstep 		= m_FParams.szPnts;
	// Get the maximum work group size for executing the kernel on the device ///////// From https://github.com/rsnemmen/OpenCL-examples/blob/e2c34f1dfefbd265cfb607c2dd6c82c799eb322a/square_array/square.c
	status = clGetKernelWorkGroupInfo(m_Kern[FUNC_FPREFIXUP], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clGetKernelWorkGroupInfo(m_Kern[FUNC_FPREFIXSUM], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clGetKernelWorkGroupInfo(m_Kern[FUNC_TALLYLISTS], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}

	// Number of total work items, calculated here after 1st image is loaded &=> know the size.
	// NB localSize must be devisor
	// NB global_work_size must be a whole number of "Preferred work group size multiple" for Nvidia.
    // 	global_work_size = ceil((float)layerstep/(float)local_work_size) * local_work_size;
    // 	local_work_size  = 32; // trial for nvidia
    // 																				if(verbosity>1){
    // 																					cout<<"\nglobal_work_size="<<global_work_size<<", local_work_size="<<local_work_size<<", deviceId="<<deviceId<<"\n"<<flush;
    // 																					cout<<"\nlayerstep=width*height="<<layerstep<<",\tsizeof(layerstep)="<< sizeof(layerstep) <<",\tsizeof(int)="<< sizeof(int) <<flush;
    // 																					cout<<"\n";
    // 																				}
	cl_int res;
	m_FParamDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(m_FParams), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	m_FluidDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(m_Fluid), 0, &res);	    if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	m_FluidTempDevice	= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(m_FluidTemp), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	m_FGenomeDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(m_FGenome), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
    m_FPrefixDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(FPrefix), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}



																				if(verbosity>1) {
																					cout << "m_FParamDevice = " 		<< m_FParamDevice << endl;
																					cout << "m_FluidDevice = " 		<< m_FluidDevice << endl;
																					cout << "m_FluidTempDevice = " 	<< m_FluidTempDevice << endl;
																					cout << "m_FGenomeDevice = " 	<< m_FGenomeDevice << endl;
																				}

	status = clEnqueueWriteBuffer(uload_queue, m_FParamDevice, 		CL_FALSE, 0, sizeof(m_FParams), 		&m_FParams, 			0, NULL, &writeEvt);	  if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.3\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, m_FluidDevice, 		CL_FALSE, 0, sizeof(m_Fluid), 		&m_Fluid, 			0, NULL, &writeEvt);	          if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.4\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, m_FluidTempDevice, 	CL_FALSE, 0, sizeof(m_FluidTemp), 				&m_FluidTemp, 		0, NULL, &writeEvt);  if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.5\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, m_FGenomeDevice, 	CL_FALSE, 0, sizeof(m_FGenome)*3, 		&m_FGenome, 0, NULL, &writeEvt);	          if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.6\n" << endl; exit_(status);}

	clFlush(uload_queue); status = clFinish(uload_queue); 					if (status != CL_SUCCESS)	{ cout << "\nclFinish(uload_queue)="				<< status << checkerror(status)<<"\n"<<flush; exit_(status);}
                                                                            if(verbosity>0) cout << "-----InitializeOpenCL finished-----\n\n" << flush;
}

void FluidSystem::CleanUp()
{
																																			cout<<"\nRunCL::CleanUp_chk0"<<flush;
	cl_int status;
	status = clReleaseMemObject(m_FParamDevice);	if (status != CL_SUCCESS)	{ cout << "\nbasemem  status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.1"<<flush;
	status = clReleaseMemObject(m_FluidDevice);	    if (status != CL_SUCCESS)	{ cout << "\nimgmem   status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.2"<<flush;
	status = clReleaseMemObject(m_FluidTempDevice);	if (status != CL_SUCCESS)	{ cout << "\ncdatabuf status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.3"<<flush;
	status = clReleaseMemObject(m_FGenomeDevice);	if (status != CL_SUCCESS)	{ cout << "\nhdatabuf status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.4"<<flush;
																																			cout<<"\nRunCL::CleanUp_chk1_finished"<<flush;
}

void FluidSystem::exit_(cl_int res)
{
	CleanUp();
	//~RunCL(); Never need to call a destructor manually.
	exit(res);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void FluidSystem::Initialize(){             //Left aside for now, implement by copying from InitializeOpenCLused for CPU only for "check_demo".

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

    size_t global_work_size = m_FParams.numBlocks * m_FParams.numThreads;
    size_t local_work_size = m_FParams.numThreads;
    cout << "-----Initialization successful-----" << endl;

    if (m_FParams.debug>1)std::cout << "Chk1.4 \n";

    // Allocate the sim parameters CUCLCUCL
    AllocateBuffer ( FPARAMS,		sizeof(FParams),	0,	1,	 GPU_OFF,     CPU_YES );//AllocateBuffer ( int buf_id, int stride,     int cpucnt, int gpucnt,    int gpumode,    int cpumode )

    cout << "----- Buffers allocated -----" << endl;

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
            clEnqueueWriteBuffer(m_queue, dest_buf, CL_TRUE, 0, cpucnt*stride, src_buf, 0, NULL, NULL);
            free(src_buf);
        }
        setBuf(&m_Fluid, buf_id, dest_buf); // stores pointer to buffer in mcpu[buf_id]
}*/

void FluidSystem::UpdateGenome (){              // Update Genome on GPU
    clCheck ( clEnqueueWriteBuffer(m_queue, m_FGenomeDevice, CL_TRUE, 0, sizeof(FGenome), &m_FGenome, 0, NULL, NULL), "FluidUpdateGenome", "clEnqueueWriteBuffer", "m_FGenomeDevice", mbDebug);
}

FGenome	FluidSystem::GetGenome(){
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

void FluidSystem::FluidParamCL (float ss, float sr, float pr, float mass, float rest, float3 bmin, float3 bmax, float estiff, float istiff, float visc, float surface_tension, float damp, float fmin, float fmax, float ffreq, float gslope, float gx, float gy, float gz, float al, float vl, float a_f, float a_p ){
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
    clCheck ( clEnqueueWriteBuffer(m_queue, m_FParamDevice, CL_TRUE, 0, sizeof(FParams), &m_FParams, 0, NULL, NULL), "FluidParamCL", "clEnqueueWriteBuffer", "m_FParamDevice", mbDebug);

}


void FluidSystem::UpdateParams (){
    // Update Params on GPU
    Vector3DF grav = Vector3DF_multiplyFloat(&m_Vec[PPLANE_GRAV_DIR], m_Param[PGRAV]);
    /*FluidParamCL (runcl,  m_Param[PSIMSCALE], m_Param[PSMOOTHRADIUS], m_Param[PRADIUS], m_Param[PMASS], m_Param[PRESTDENSITY],
                      *(float3*)& m_Vec[PBOUNDMIN], *(float3*)& m_Vec[PBOUNDMAX], m_Param[PEXTSTIFF], m_Param[PINTSTIFF],
                      m_Param[PVISC], m_Param[PSURFACE_TENSION], m_Param[PEXTDAMP], m_Param[PFORCE_MIN], m_Param[PFORCE_MAX], m_Param[PFORCE_FREQ],
                      m_Param[PGROUND_SLOPE], grav.x, grav.y, grav.z, m_Param[PACCEL_LIMIT], m_Param[PVEL_LIMIT],
                      m_Param[PACTUATION_FACTOR], m_Param[PACTUATION_PERIOD]);*/
}

void FluidSystem::SetParam (int p, float v ){
    m_Param[p] = v;
    UpdateParams ();
}

void FluidSystem::SetVec (int p, Vector3DF v ){
    m_Vec[p] = v;
    UpdateParams ();
}
/*      Probably unnecessary, no desctroctor needed in OpenCL, rest is done in RunCL::CleanUp()
void FluidSystem::Exit (){
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

void FluidSystem::AllocateBuffer ( int buf_id, int stride, int cpucnt, int gpucnt, int gpumode, int cpumode ){   // mallocs a buffer - called by FluidSystem::Initialize(), AllocateParticles, and AllocateGrid()
//also called by WriteDemoSimParams(..)
    cl_int err;
    bool rtn = true;
    if (m_FParams.debug>1)std::cout<<"\nAllocateBuffer ( int buf_id="<<buf_id<<", int stride="<<stride<<", int cpucnt="<<cpucnt<<", int gpucnt="<<gpucnt<<", int "<<gpumode<<", int "<<cpumode<<" )\t"<<std::flush;
    if (cpumode == CPU_YES) {
        char* src_buf  = bufC(&m_Fluid, buf_id);
        char* dest_buf = (char*) malloc(cpucnt*stride); //  ####  malloc the buffer   ####
        if (src_buf != 0x0) {
            memcpy(dest_buf, src_buf, cpucnt*stride);
            free(src_buf);
        }
        setBuf(&m_Fluid, buf_id, dest_buf); // stores pointer to buffer in mcpu[buf_id]
}

    if(gpumode == GPU_SINGLE || gpumode == GPU_DUAL || gpumode == GPU_TEMP){
        cl_int err;
        cl_ulong free1, free2, total;
        clGetDeviceInfo(m_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(total), &total, NULL);
        //if (m_FParams.debug>1)printf("\nOpenCL Memory: free=%lu, total=%lu.\t",free1,total);

        if (gpumode == GPU_SINGLE || gpumode == GPU_DUAL )    {
            if (gpuptr(&m_Fluid, buf_id) != 0x0)            clCheck(clReleaseMemObject(gpuVar(&m_Fluid, buf_id)), "AllocateBuffer", "clReleaseMemObject", "Fluid.gpu", mbDebug);
                                                            cl_int err;
                                                            setGpuBuf(&m_Fluid, buf_id, clCreateBuffer(m_context, CL_MEM_READ_WRITE, stride*gpucnt, NULL, &err));

            if (m_FParams.debug>1)                          std::cout<<"\t\t gpuptr(&m_Fluid, "<<buf_id<<")'"<<gpuptr(&m_Fluid, buf_id)<<",   gpuVar(&m_Fluid, "<<buf_id<<")="<<gpuVar(&m_Fluid, buf_id)<<"\t"<<std::flush;
            //if(err != CL_SUCCESS)                           FluidSystem::Exit();
        }

        if (gpumode == GPU_TEMP || gpumode == GPU_DUAL ) {
            if (gpuptr(&m_FluidTemp, buf_id) != 0x0)        clCheck(clReleaseMemObject(gpuVar(&m_FluidTemp, buf_id)), "AllocateBuffer", "clReleaseMemObject", "FluidTemp.gpu", mbDebug);

                                                            setGpuBuf(&m_FluidTemp, buf_id, clCreateBuffer(m_context, CL_MEM_READ_WRITE, stride*gpucnt, NULL, &err));
            //if (err != CL_SUCCESS)                          FluidSystem::Exit();
        }

        //no clFinish needed, as EnqueueWriteBuffer function gets called with blocking_wirte set to CL_TRUE
        if (m_FParams.debug>1)printf("\nAfter allocation: free=%lu, total=%lu, this buffer=%lu.\n",free2,total,(free1-free2) );
    }
}

void FluidSystem::AllocateParticles ( int cnt, int gpu_mode, int cpu_mode ){ // calls AllocateBuffer(..) for each buffer.
// Defaults in header : int gpu_mode = GPU_DUAL, int cpu_mode = CPU_YES
// Called by FluidSystem::ReadPointsCSV(..), and FluidSystem::WriteDemoSimParams(...), cnt = mMaxPoints.
    size_t threefloat = 3 * sizeof(float);
if (m_FParams.debug>1)std::cout<<"\n\nAllocateParticles ( int cnt="<<cnt<<", int "<<gpu_mode<<", int "<<cpu_mode<<" ), debug="<<m_FParams.debug<<", launchParams.debug="<<launchParams.debug<<"\t";//<<std::flush;
if (m_FParams.debug>1)std::cout<<"\tGPU_OFF=0, GPU_SINGLE=1, GPU_TEMP=2, GPU_DUAL=3, CPU_OFF=4, CPU_YES=5"<<std::flush;
    AllocateBuffer ( FPOS,		sizeof(Vector3DF),	cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FCLR,		sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FVEL,		sizeof(Vector3DF),	cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FVEVAL,	sizeof(Vector3DF),	cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FAGE,		sizeof(uint),       cnt,    m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FPRESS,	sizeof(float),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FDENSITY,	sizeof(float),		cnt, 	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FFORCE,	sizeof(threefloat),	cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FCLUSTER,	sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FGCELL,	sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FGNDX,		sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FGNEXT,	sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FNBRNDX,	sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FNBRCNT,	sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FSTATE,	sizeof(uint),		cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    // extra buffers for morphogenesis
    AllocateBuffer ( FELASTIDX,	    sizeof(uint[BOND_DATA]),             cnt,   m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FPARTICLEIDX,	sizeof(uint[BONDS_PER_PARTICLE *2]), cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FPARTICLE_ID,	sizeof(uint),		                 cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FMASS_RADIUS,	sizeof(uint),		                 cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FNERVEIDX,	    sizeof(uint),		                 cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FCONC,	        sizeof(float[NUM_TF]),		         cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FEPIGEN,	    sizeof(uint[NUM_GENES]),	         cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    //AllocateBuffer ( FCURAND_STATE,	sizeof(curandState_t),	             cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );
    AllocateBuffer ( FCURAND_SEED,	sizeof(unsigned long long),	         cnt,	m_FParams.szPnts,	gpu_mode, cpu_mode );

    // Update GPU access pointers
    if (gpu_mode != GPU_OFF ) {

        clCheck( clEnqueueWriteBuffer(m_queue, m_FluidDevice, CL_TRUE, 0, sizeof(FBufs), &m_Fluid, 0, NULL, NULL), "AllocateParticles", "clEnqueueWriteBuffer", "m_FluidDevice", mbDebug);
        clCheck( clEnqueueWriteBuffer(m_queue, m_FluidTempDevice, CL_TRUE, 0, sizeof(FBufs), &m_FluidTemp, 0, NULL, NULL),	"AllocateParticles", "clEnqueueWriteBuffer", "m_FluidTempDevice", mbDebug);
        clCheck( clEnqueueWriteBuffer(m_queue, m_FParamDevice, CL_TRUE, 0, sizeof(FParams), &m_FParams, 0, NULL, NULL),  "AllocateParticles", "clEnqueueWriteBuffer", "m_FParamDevice", mbDebug);
        clCheck( clEnqueueWriteBuffer(m_queue, m_FGenomeDevice, CL_TRUE, 0, sizeof(FGenome), &m_FGenome, 0, NULL, NULL),  "AllocateParticles", "clEnqueueWriteBuffer", "m_FGenomeDevice", mbDebug);

        clCheck(clFinish(m_queue), "AllocateParticles", "clFinish", "", mbDebug );
    }

    // Allocate auxiliary buffers (prefix sums)
    int blockSize = SCAN_BLOCKSIZE << 1;
    int numElem1 = m_GridTotal;
    int numElem2 = int ( numElem1 / blockSize ) + 1;
    int numElem3 = int ( numElem2 / blockSize ) + 1;

    if (gpu_mode != GPU_OFF ) {
        AllocateBuffer ( FAUXARRAY1,	sizeof(uint),		0,	numElem2, GPU_SINGLE, CPU_OFF );
        AllocateBuffer ( FAUXSCAN1,	    sizeof(uint),		0,	numElem2, GPU_SINGLE, CPU_OFF );
        AllocateBuffer ( FAUXARRAY2,	sizeof(uint),		0,	numElem3, GPU_SINGLE, CPU_OFF );
        AllocateBuffer ( FAUXSCAN2,	    sizeof(uint),		0,	numElem3, GPU_SINGLE, CPU_OFF );
    }
}

void FluidSystem::AllocateBufferDenseLists(int buf_id, int stride, int gpucnt, int lists) {
    // mallocs a buffer - called by FluidSystem::AllocateGrid(int gpu_mode, int cpu_mode)
    // Need to save "pointers to the allocated gpu buffers" in a cpu array, AND then clEnqueueWriteBuffer(...) that list of pointers into the device array.
    // also called by FluidSystem::....() to quadruple buffer as needed.
    clCheck(clFinish(m_queue), "AllocateBufferDenseLists ", "clFinish", "before 1st clGetDeviceInfo", mbDebug);
    cl_ulong total;
    clGetDeviceInfo(m_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &total, NULL);
    if (m_FParams.debug > 1) printf("\nOpenCL Memory: total=%lu.\t", total);

    cl_mem* listpointer = (cl_mem*)&bufC(&m_Fluid, lists)[buf_id * sizeof(cl_mem)];

    if (m_FParams.debug > 1) printf("\n*listpointer=%p, listpointer=%p,  lists=%i, buf_id=%i, \t", (cl_mem*)*listpointer, listpointer, lists, buf_id);

    if (m_FParams.debug > 1) printf("\nAllocateBufferDenseLists: buf_id=%i, stride=%i, gpucnt=%i, lists=%i,  .\t", buf_id, stride, gpucnt, lists);

    if (*listpointer != NULL) clCheck(clReleaseMemObject(*listpointer), "AllocateBufferDenseLists1", "clReleaseMemObject", "*listpointer", mbDebug);

    cl_int err;
    *listpointer = clCreateBuffer(m_context, CL_MEM_READ_WRITE, stride * gpucnt, NULL, &err);
    bool result = clCheck(err, "AllocateBufferDenseLists2", "clCreateBuffer", "listpointer", mbDebug);

    clCheck(clFinish(m_queue), "AllocateBufferDenseLists ", "clFinish", "after allocation", mbDebug);
    if (result == false) _exit(result);
}


// void FluidSystem::AllocateBufferDenseLists ( int buf_id, int stride, int gpucnt, int lists ) {    // mallocs a buffer - called by FluidSystem::AllocateGrid(int gpu_mode, int cpu_mode)
// // Need to save "pointers to the allocated gpu buffers" in a cpu array, AND then cuMemcpyHtoD(...) that list of pointers into the device array.
//     // also called by FluidSystem::....()  to quadruple buffer as needed.
//     clCheck(clFinish(m_queue), "AllocateBufferDenseLists ", "clFinish", "before 1st cudaMemGetInfo(&free1, &total)", mbDebug);
//     size_t   free1, free2, total;
//     cudaMemGetInfo(&free1, &total);
//     if (m_FParams.debug>1)printf("\nCuda Memory: free=%lu, total=%lu.\t",free1,total);
//
//     cl_device_idptr*  listpointer = (cl_device_idptr*) &bufC(&m_Fluid, lists)[buf_id * sizeof(cl_device_idptr)] ;
//
//     //cl_device_idptr  listpointer2 = m_Fluid.gpuptr(lists)[buf_id]  ;
//     if (m_FParams.debug>1)printf("\n*listpointer=%p, listpointer=%p,  lists=%i, buf_id=%i, \t", (cl_device_idptr* ) *listpointer, listpointer,  /*listpointer2,*/ lists, buf_id);/*listpointer2=%llu,*/
//     //if (m_FParams.debug>1)cout<<"\n listpointer is an:"<< typeid(listpointer).name()<<" *listpointer is an:"<< typeid(*listpointer).name()<<" listpointer2 is an:"<< typeid(listpointer2).name()<<" .  "<<std::flush;//" *listpointer2 is an:"<< typeid(*listpointer2).name()<<
//
//     if (m_FParams.debug>1)printf("\nAllocateBufferDenseLists: buf_id=%i, stride=%i, gpucnt=%i, lists=%i,  .\t", buf_id, stride, gpucnt, lists);
//     if (*listpointer != 0x0) clCheck(cuMemFree(*listpointer), "AllocateBufferDenseLists1", "cuMemFree", "*listpointer", mbDebug);
//     bool result = clCheck( cuMemAlloc( listpointer, stride*gpucnt),   "AllocateBufferDenseLists2", "cuMemAlloc", "listpointer", mbDebug);
//
//     clCheck(clFinish(m_queue), "AllocateBufferDenseLists ", "clFinish", "before 2nd cudaMemGetInfo(&free2, &total)", mbDebug);
//     cudaMemGetInfo(&free2, &total);
//     if (m_FParams.debug>1)printf("\nAfter allocation: free=%lu, total=%lu, this buffer=%lu.\n",free2,total,(free1-free2) );
//     if(result==false)Exit();
// }

// void FluidSystem::AllocateGrid(int gpu_mode, int cpu_mode){ // NB void FluidSystem::AllocateBuffer (int buf_id, int stride, int cpucnt, int gpucnt, int gpumode, int cpumode)
//     // Allocate grid
//     int cnt = m_GridTotal;
//     m_FParams.szGrid = (m_FParams.gridBlocks * m_FParams.gridThreads);
//     if (m_FParams.debug>1)cout<<"\nAllocateGrid: m_FParams.szGrid = ("<<m_FParams.gridBlocks<<" * "<<m_FParams.gridThreads<<")"<<std::flush;
//     AllocateBuffer ( FGRID,		sizeof(uint),		mMaxPoints,	m_FParams.szPnts,	gpu_mode, cpu_mode );    // # grid elements = number of points
//     AllocateBuffer ( FGRIDCNT,	sizeof(uint),		cnt,	    m_FParams.szGrid,	gpu_mode, cpu_mode );
//     AllocateBuffer ( FGRIDOFF,	sizeof(uint),		cnt,	    m_FParams.szGrid,	gpu_mode, cpu_mode );
//     AllocateBuffer ( FGRIDACT,	sizeof(uint),		cnt,	    m_FParams.szGrid,	gpu_mode, cpu_mode );      // ?? not used ?? ... active bins i.e. containing particles ?
//     // extra buffers for dense lists
//     AllocateBuffer ( FGRIDCNT_ACTIVE_GENES,  sizeof(uint[NUM_GENES]),       cnt,   m_FParams.szGrid,	gpu_mode, cpu_mode );
//     AllocateBuffer ( FGRIDOFF_ACTIVE_GENES,  sizeof(uint[NUM_GENES]),       cnt,   m_FParams.szGrid,	gpu_mode, cpu_mode );
//     AllocateBuffer ( FDENSE_LIST_LENGTHS,	 sizeof(uint),		      NUM_GENES,   NUM_GENES,	        gpu_mode, cpu_mode );
//     AllocateBuffer ( FDENSE_LISTS,	         sizeof(cl_mem),     NUM_GENES,   NUM_GENES,           gpu_mode, cpu_mode );
//     AllocateBuffer ( FDENSE_BUF_LENGTHS,	 sizeof(uint),            NUM_GENES,   NUM_GENES,           gpu_mode, cpu_mode );
//
//     AllocateBuffer ( FGRIDCNT_CHANGES,               sizeof(uint[NUM_CHANGES]),       cnt,   m_FParams.szGrid,	    gpu_mode, cpu_mode );
//     AllocateBuffer ( FGRIDOFF_CHANGES,               sizeof(uint[NUM_CHANGES]),       cnt,   m_FParams.szGrid,	    gpu_mode, cpu_mode );
//     AllocateBuffer ( FDENSE_LIST_LENGTHS_CHANGES,	 sizeof(uint),		      NUM_CHANGES,   NUM_CHANGES,	        gpu_mode, cpu_mode );
//     AllocateBuffer ( FDENSE_LISTS_CHANGES,	         sizeof(cl_mem),     NUM_CHANGES,   NUM_CHANGES,           gpu_mode, cpu_mode );
//     AllocateBuffer ( FDENSE_BUF_LENGTHS_CHANGES,	 sizeof(uint),            NUM_CHANGES,   NUM_CHANGES,           gpu_mode, cpu_mode );
//
//     if (gpu_mode != GPU_OFF ) {
//         /*if(gpu_mode == GPU_SINGLE || gpu_mode == GPU_DUAL )*/
//         cl_int status;
//
//         for(int i=0; i<NUM_GENES; i++){ //for each gene allocate intial buffer, write pointer and size to FDENSE_LISTS and FDENSE_LIST_LENGTHS
//             cl_device_idptr*  _listpointer = (cl_device_idptr*) &bufC(&m_Fluid, FDENSE_LISTS)[i * sizeof(cl_device_idptr)] ;
//             *_listpointer = 0x0;
//             AllocateBufferDenseLists( i, sizeof(uint), INITIAL_BUFFSIZE_ACTIVE_GENES, FDENSE_LISTS);  // AllocateBuffer writes pointer to  m_Fluid.gpuptr(buf_id).
//             bufI(&m_Fluid, FDENSE_LIST_LENGTHS)[i] = 0;
//             bufI(&m_Fluid, FDENSE_BUF_LENGTHS)[i]  = INITIAL_BUFFSIZE_ACTIVE_GENES;
//         }
//         /*if(gpu_mode == GPU_SINGLE || gpu_mode == GPU_DUAL )*/
//         for(int i=0; i<NUM_CHANGES; i++){ //Same for the changes lists
//             cl_device_idptr*  _listpointer = (cl_device_idptr*) &bufC(&m_Fluid, FDENSE_LISTS_CHANGES)[i * sizeof(cl_device_idptr)] ;
//             *_listpointer = 0x0;
//             AllocateBufferDenseLists( i, sizeof(uint), 2*INITIAL_BUFFSIZE_ACTIVE_GENES, FDENSE_LISTS_CHANGES); // NB buf[2][list_length] holding : particleIdx, bondIdx
//             bufI(&m_Fluid, FDENSE_LIST_LENGTHS_CHANGES)[i] = 0;
//             bufI(&m_Fluid, FDENSE_BUF_LENGTHS_CHANGES)[i]  = INITIAL_BUFFSIZE_ACTIVE_GENES;
//         }
//
//         clEnqueueWriteBuffer(m_queue, gpuVar(&m_Fluid, FDENSE_LISTS), CL_TRUE, 0, NUM_GENES * sizeof(cl_mem), bufC(&m_Fluid, FDENSE_LISTS), 0, NULL, NULL);
//         clEnqueueWriteBuffer(m_queue, gpuVar(&m_Fluid, FDENSE_LISTS_CHANGES), CL_TRUE, 0, NUM_GENES * sizeof(cl_mem), bufC(&m_Fluid, FDENSE_LISTS_CHANGES), 0, NULL, NULL);
//         clEnqueueWriteBuffer(m_queue, gpuVar(&m_Fluid, FDENSE_BUF_LENGTHS), CL_TRUE, 0, NUM_GENES * sizeof(cl_mem), bufC(&m_Fluid, FDENSE_BUF_LENGTHS), 0, NULL, NULL);
//
//         // Copy FDENSE_BUF_LENGTHS_CHANGES
//         status = clEnqueueWriteBuffer(m_queue, gpuVar(&m_Fluid, FDENSE_BUF_LENGTHS_CHANGES), CL_TRUE, 0, sizeof(uint[NUM_CHANGES]), bufI(&m_Fluid, FDENSE_BUF_LENGTHS_CHANGES), 0, NULL, NULL);
//         if (status != CL_SUCCESS) {if (mbDebug) {printf("Error in AllocateGrid: clEnqueueWriteBuffer for FDENSE_BUF_LENGTHS_CHANGES failed with status %d\n", status);}}
//
//         // Copy m_FluidDevice
//         status = clEnqueueWriteBuffer(m_queue, m_FluidDevice, CL_TRUE, 0, sizeof(FBufs), &m_Fluid, 0, NULL, NULL);// Update GPU access pointers
//         if (status != CL_SUCCESS) {if (mbDebug) {printf("Error in AllocateGrid: clEnqueueWriteBuffer for m_FluidDevice failed with status %d\n", status);}}
//
//         clCheck(clFinish(m_queue), "AllocateParticles", "clFinish", " ", mbDebug);
//     }
// }

void FluidSystem::SavePointsCSV2 ( const char * relativePath, int frame ){
    if (m_FParams.debug>1) std::cout << "\n  SavePointsCSV2 ( const char * relativePath = "<< relativePath << ", int frame = "<< frame << " );  started \n" << std::flush;
    char buf[256];
    frame += 100000;    // ensures numerical and alphabetic order match
    sprintf ( buf, "%s/particles_pos_vel_color%04d.csv", relativePath, frame );
    FILE* fp = fopen ( buf, "w" );
    if (fp == NULL) {
        if (m_FParams.debug>1) std::cout << "\nvoid FluidSystem::SavePointsCSV ( const char * relativePath, int frame )  Could not open file "<< fp <<"\n"<< std::flush;
        assert(0);
    }
    int numpnt = mMaxPoints;//NumPoints();
    Vector3DF* Pos;
    Vector3DF* Vel;
    float *Conc;
    uint* Age, *Clr, *NerveIdx, *ElastIdx, *Particle_Idx, *Particle_ID, *Mass_Radius, *EpiGen;                  // Q: why are these pointers? A: they get dereferenced below.
    uint mass, radius;
    float *ElastIdxPtr;

    fprintf(fp, "i,, x coord, y coord, z coord\t\t x vel, y vel, z vel\t\t age,  color\t\t FELASTIDX[%u*%u]", BONDS_PER_PARTICLE, DATA_PER_BOND);  // This system inserts commas to align header with csv data
    for (int i=0; i<BONDS_PER_PARTICLE; i++)fprintf(fp, ",(%u)[0]curIdx, [1]elastLim, [2]restLn, [3]modulus, [4]damping, [5]partID, [6]stress_sq integrator, [7]stress integrator, [8]change-type,,  ",i);
    fprintf(fp, "\t");
    fprintf(fp, "\tParticle_ID, mass, radius, FNERVEIDX,\t\t Particle_Idx[%u*2]", BONDS_PER_PARTICLE);
    for (int i=0; i<BONDS_PER_PARTICLE; i++)fprintf(fp, "%u,,, ",i);
    fprintf(fp, "\t\tFCONC[%u] ", NUM_TF);
    for (int i=0; i<NUM_TF; i++)fprintf(fp, "%u, ",i);
    fprintf(fp, "\t\tFEPIGEN[%u] ", NUM_GENES);
    for (int i=0; i<NUM_GENES; i++)fprintf(fp, "%u, ",i);
    fprintf(fp, "\n");

    for(int i=0; i<numpnt; i++) {       // nb need get..() accessors for private data.
        Pos = getPos(i);                // e.g.  Vector3DF* getPos ( int n )	{ return &m_Fluid.bufV3(FPOS)[n]; }
        Vel = getVel(i);
        Age = getAge(i);
        Clr = getClr(i);
        ElastIdx = getElastIdx(i);      // NB [BONDS_PER_PARTICLE]
      //if (m_FParams.debug>1)printf("\t%u,",ElastIdx[0]);
        ElastIdxPtr = (float*)ElastIdx; // #############packing floats and uints into the same array - should replace with a struct.#################
        Particle_Idx = getParticle_Idx(i);
        Particle_ID = getParticle_ID(i);//# uint  original pnum, used for bonds between particles. 32bit, track upto 4Bn particles.
        if(*Particle_ID==0){
         if (m_FParams.debug>1) std::cout << "SavePointsCSV2: Particle_ID = pointer not assigned. i="<<i<<". \t" << std::flush;
         return;
        }
        // ? should I be splitting mass_radius with bitshift etc  OR just use two uit arrays .... where are/will these used anyway ?
        Mass_Radius = getMass_Radius(i);//# uint holding modulus 16bit and limit 16bit.
        if(*Mass_Radius==0){   mass = 0; }else{  mass = *Mass_Radius; }    // modulus          // '&' bitwise AND is bit masking. ;
        radius = mass >> 16;
        mass = mass & TWO_POW_16_MINUS_1;

        NerveIdx = getNerveIdx(i);      //# uint
        //Conc = getConc(i);              //# float[NUM_TF]        NUM_TF = num transcription factors & morphogens
        //EpiGen = getEpiGen(i);          //# uint[NUM_GENES]  see below.

        fprintf(fp, "%u,,%f,%f,%f,\t%f,%f,%f,\t %u, %u,, \t", i, Pos->x, Pos->y,Pos->z, Vel->x,Vel->y,Vel->z, *Age, *Clr );
        //if (m_FParams.debug>1) std::cout<<"\t"<<Pos->z<<std::flush;
        for(int j=0; j<(BOND_DATA); j+=DATA_PER_BOND) {
            fprintf(fp, "%u, %f, %f, %f, %f, %u, %f, %f, %u, ", ElastIdx[j], ElastIdxPtr[j+1], ElastIdxPtr[j+2], ElastIdxPtr[j+3], ElastIdxPtr[j+4], ElastIdx[j+5], ElastIdxPtr[j+6], ElastIdxPtr[j+7], ElastIdx[j+8] );

           /*
            // if ((j%DATA_PER_BOND==0)||((j+1)%DATA_PER_BOND==0))  fprintf(fp, "%u, ",  ElastIdx[j] );  // print as int   [0]current index, [5]particle ID, [6]bond index
           // else  fprintf(fp, "%f, ",  ElastIdxPtr[j] );                                               // print as float [1]elastic limit, [2]restlength, [3]modulus, [4]damping coeff,
           //  if((j+1)%DATA_PER_BOND==0)
            */
            fprintf(fp, "\t\t");
        }
        fprintf(fp, " \t%u, %u, %u, %u, \t\t", *Particle_ID, mass, radius, *NerveIdx );
        for(int j=0; j<(BONDS_PER_PARTICLE*2); j+=2)   { fprintf(fp, "%u, %u,, ",  Particle_Idx[j], Particle_Idx[j+1] );}  fprintf(fp, "\t\t"); // NB index of other particle AND other particle's index of the bond

        for(int j=0; j<(NUM_TF); j++)               {
            Conc = getConc(j);
            fprintf(fp, "%f, ",  Conc[i] );
        }fprintf(fp, "\t\t");

        for(int j=0; j<(NUM_GENES); j++)            {
            EpiGen = getEpiGen(j);
            fprintf(fp, "%u, ",  EpiGen[i] );   // NB FEPIGEN[gene][particle], for memory efficiency on the device. ? Need to test.
        }fprintf(fp, " \n");
    }
    fclose ( fp );
    fflush ( fp );
}

void FluidSystem::Run2PhysicalSort(){
                                                                                if(verbosity>0) cout << "-----starting Run2PhysicalSort-----\n\n" << flush;

    if(m_FParams.debug>1)std::cout<<"\n####\nRun2PhysicalSort()start";
    InsertParticlesCL ( 0x0, 0x0, 0x0 );
    clCheck(clFinish(m_queue), "Run", "clFinish", "After InsertParticlesCL", mbDebug);
    if(launchParams.debug>0){
        std::cout<<"\nchk a"<<std::flush;
        TransferFromCL ();
        m_Debug_file++;
        std::cout<<"\nchk b: launchParams.outPath="<<launchParams.outPath<<",  m_Frame+m_Debug_file=0;="<<m_Frame+m_Debug_file<<"\t"<<std::flush;
        SavePointsCSV2 (  launchParams.outPath, m_Frame+m_Debug_file );
        std::cout << "\n\nRun2PhysicalSort() Chk1, saved "<<launchParams.outPath<< m_Frame+m_Debug_file <<".csv  After  InsertParticlesCL\n"<<std::flush;
        //TransferFromTempCL(int buf_id, int sz );
    }
    PrefixSumCellsCL ( 1 );
    clCheck(clFinish(m_queue), "Run", "clFinish", "After PrefixSumCellsCL", mbDebug);
    if(launchParams.debug>0){
        TransferFromCL ();
        m_Debug_file++;
        SavePointsCSV2 (  launchParams.outPath, m_Frame+m_Debug_file );
        std::cout << "\n\nRun2PhysicalSort() Chk2, saved "<<launchParams.outPath<< m_Frame+m_Debug_file <<".csv  After  PrefixSumCellsCL\n"<<std::flush;
        //TransferFromTempCL(int buf_id, int sz );
    }
    CountingSortFullCL ( 0x0 );
    clCheck(clFinish(m_queue), "Run", "clFinish", "After CountingSortFullCL", mbDebug);
    if(m_FParams.debug>1)std::cout<<"\n####\nRun2PhysicalSort()end";
}

void FluidSystem::AdvanceTime () {  // may need to prune unused details from this fn.
    m_Time += m_DT;
}

///////////////////////////////////////////////////
// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale)
void FluidSystem::SetupGrid ( Vector3DF min, Vector3DF max, float sim_scale, float cell_size){
    float world_cellsize = cell_size / sim_scale;
    m_GridMin = min;
    m_GridMax = max;
    m_GridSize = m_GridMax;
    m_GridSize = Vector3DF_subtractVector3DF(&m_GridSize,&m_GridMin);
    //m_GridSize -= m_GridMin;
    m_GridRes.x = (int) ceil ( m_GridSize.x / world_cellsize );		// Determine grid resolution
    m_GridRes.y = (int) ceil ( m_GridSize.y / world_cellsize );
    m_GridRes.z = (int) ceil ( m_GridSize.z / world_cellsize );
    m_GridSize.x = m_GridRes.x * cell_size / sim_scale;				// Adjust grid size to multiple of cell size
    m_GridSize.y = m_GridRes.y * cell_size / sim_scale;
    m_GridSize.z = m_GridRes.z * cell_size / sim_scale;
    m_GridDelta = *Vector3DF_operator_equal_Vector3DI(&m_GridDelta, &m_GridRes);		// delta = translate from world space to cell #
    m_GridDelta = Vector3DF_divideByVector3DF(&m_GridDelta, m_GridSize);
    m_GridTotal = (int)(m_GridRes.x * m_GridRes.y * m_GridRes.z);

    // Number of cells to search:
    // n = (2r / w) +1,  where n = 1D cell search count, r = search radius, w = world cell width
    m_GridSrch = (int) (floor(2.0f*(m_Param[PSMOOTHRADIUS]/sim_scale) / world_cellsize) + 1.0f);
    if ( m_GridSrch < 2 ) m_GridSrch = 2;
    m_GridAdjCnt = m_GridSrch * m_GridSrch * m_GridSrch ;			// 3D search count = n^3, e.g. 2x2x2=8, 3x3x3=27, 4x4x4=64

    if ( m_GridSrch > 6 ) {
        //if (m_FParams.debug>1)nvprintf ( "ERROR: Neighbor search is n > 6. \n " );
        exit(-1);
    }

    int cell = 0;
    for (int y=0; y < m_GridSrch; y++ )
        for (int z=0; z < m_GridSrch; z++ )
            for (int x=0; x < m_GridSrch; x++ )
                m_GridAdj[cell++] = ( y*m_GridRes.z + z )*m_GridRes.x +  x ;			// -1 compensates for ndx 0=empty

    if ( mPackGrid != 0x0 ) free ( mPackGrid );
    mPackGrid = (int*) malloc ( sizeof(int) * m_GridTotal );
}

void FluidSystem::SetupSPH_Kernels (){
    m_Param [ PDIST ] = pow ( (float) m_Param[PMASS] / m_Param[PRESTDENSITY], 1.0f/3.0f );
    m_R2 = m_Param [PSMOOTHRADIUS] * m_Param[PSMOOTHRADIUS];
    m_Poly6Kern = 315.0f / (64.0f * 3.141592f * pow( m_Param[PSMOOTHRADIUS], 9.0f) );	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
    m_SpikyKern = -45.0f / (3.141592f * pow( m_Param[PSMOOTHRADIUS], 6.0f) );			// Laplacian of viscocity (denominator): PI h^6
    m_LapKern = 45.0f / (3.141592f * pow( m_Param[PSMOOTHRADIUS], 6.0f) );
}

void FluidSystem::SetupDefaultParams (){
    //  Range = +/- 10.0 * 0.006 (r) =	   0.12			m (= 120 mm = 4.7 inch)
    //  Container Volume (Vc) =			   0.001728		m^3
    //  Rest Density (D) =				1000.0			kg / m^3
    //  Particle Mass (Pm) =			   0.00020543	kg						(mass = vol * density)
    //  Number of Particles (N) =		4000.0
    //  Water Mass (M) =				   0.821		kg (= 821 grams)
    //  Water Volume (V) =				   0.000821     m^3 (= 3.4 cups, .21 gals)
    //  Smoothing Radius (R) =             0.02			m (= 20 mm = ~3/4 inch)
    //  Particle Radius (Pr) =			   0.00366		m (= 4 mm  = ~1/8 inch)
    //  Particle Volume (Pv) =			   2.054e-7		m^3	(= .268 milliliters)
    //  Rest Distance (Pd) =			   0.0059		m
    //
    //  Given: D, Pm, N
    //    Pv = Pm / D			0.00020543 kg / 1000 kg/m^3 = 2.054e-7 m^3
    //    Pv = 4/3*pi*Pr^3    cuberoot( 2.054e-7 m^3 * 3/(4pi) ) = 0.00366 m
    //     M = Pm * N			0.00020543 kg * 4000.0 = 0.821 kg
    //     V =  M / D              0.821 kg / 1000 kg/m^3 = 0.000821 m^3
    //     V = Pv * N			 2.054e-7 m^3 * 4000 = 0.000821 m^3
    //    Pd = cuberoot(Pm/D)    cuberoot(0.00020543/1000) = 0.0059 m
    //
    // Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
    // Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
    //    (k = number of cells, gs = cell size, d = simulation scale)

    // "The viscosity coefficient is the dynamic viscosity, visc > 0 (units Pa.s),
    // and to include a reasonable damping contribution, it should be chosen
    // to be approximately a factor larger than any physical correct viscosity
    // coefficient that can be looked up in the literature. However, care should
    // be taken not to exaggerate the viscosity coefficient for fluid materials.
    // If the contribution of the viscosity force density is too large, the net effect
    // of the viscosity term will introduce energy into the system, rather than
    // draining the system from energy as intended."
    //    Actual visocity of water = 0.001 Pa.s    // viscosity of water at 20 deg C.

    m_Time = 0.0f;							// Start at T=0
    m_DT = 0.003f;

    m_Param [ PSIMSCALE ] =		0.005f;			// unit size
    m_Param [ PVISC ] =			0.50f;			// pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)
    m_Param [ PSURFACE_TENSION ] = 0.1f;
    m_Param [ PRESTDENSITY ] =	400.0f;			// kg / m^3
    m_Param [ PSPACING ]	=	0.0f;			// spacing will be computed automatically from density in most examples (set to 0 for autocompute)
    m_Param [ PMASS ] =			0.00020543f;		// kg
    m_Param [ PRADIUS ] =		0.015f;			// m
    m_Param [ PDIST ] =			0.0059f;			// m
    m_Param [ PSMOOTHRADIUS ] =	0.015f;			// m
    m_Param [ PINTSTIFF ] =		1.0f;
    m_Param [ PEXTSTIFF ] =		50000.0f;
    m_Param [ PEXTDAMP ] =		100.0f;
    m_Param [ PACCEL_LIMIT ] =	150.0f;			// m / s^2
    m_Param [ PVEL_LIMIT ] =	3.0f;			// m / s
    m_Param [ PGRAV ] =			1.0f;

    m_Param [ PGROUND_SLOPE ] = 0.0f;
    m_Param [ PFORCE_MIN ] =	0.0f;
    m_Param [ PFORCE_MAX ] =	0.0f;
    m_Param [ PFORCE_FREQ ] =	16.0f;
    m_Vec [ PPLANE_GRAV_DIR ] = Vector3DF_init_with_values( 0, -9.8f, 0 );

    // Default sim config
    m_Param [PGRIDSIZE] = m_Param[PSMOOTHRADIUS] * 2;

    m_Param [ PACTUATION_FACTOR ] = 0;
    m_Param [ PACTUATION_PERIOD ] = 1;
}

void FluidSystem::SetupExampleParams (uint spacing){
    Vector3DF pos;
    Vector3DF min, max;
    m_Param [ PSPACING ] = spacing;

    //std::cout<<"\nSetupExampleParams()1: m_Param[PEXAMPLE] = "<<m_Param[PEXAMPLE]<<"\n"<<std::flush;
    //std::cout<<"\nSetupExampleParams()2: launchParams.genomePath = "<<launchParams.genomePath<<"\n"<<std::flush;

    switch ( (int) m_Param[PEXAMPLE] ) {

    case 0:	{	// Regression test. N x N x N static grid

        int k = (int) ceil ( pow ( (float) m_Param[PNUM], (float) 1.0f/3.0f ) );
        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values( 0, 0, 0 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values( 2.0f+(k/2), 2.0f+(k/2), 2.0f+(k/2) );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( 1.0f, 1.0f, 1.0f );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 1.0f+(k/2), 1.0f+(k/2), 1.0f+(k/2) );

        m_Param [ PGRAV ] = 0.0;
        m_Vec [ PPLANE_GRAV_DIR ] = Vector3DF_init_with_values( 0.0, 0.0, 0.0 );
        //m_Param [ PSPACING ] = spacing;//0.5;				// Fixed spacing		Dx = x-axis density
        m_Param [ PSMOOTHRADIUS ] =	m_Param [PSPACING];		// Search radius
        //m_Toggle [ PRUN ] = false;				// Do NOT run sim. Neighbors only.
        //m_Param [PDRAWMODE] = 1;				// Point drawing
        //m_Param [PDRAWGRID] = 1;				// Grid drawing
        //m_Param [PDRAWTEXT] = 1;				// Text drawing
        m_Param [PSIMSCALE ] = 1.0f;
        launchParams.read_genome = 'y';
    }
    break;
    case 1:		// Tower
        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values(   0,   0,   0 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values(  256, 128, 256 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values(  5,   5,  5 );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 256*0.3, 128*0.9, 256*0.3 );
        break;
    case 2:		// Wave pool
        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values(   0,   0,   0 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values(  400, 200, 400 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( 100, 80,  100 );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 300, 190, 300 );
        m_Param [ PFORCE_MIN ] = 100.0f;
        m_Param [ PFORCE_FREQ ] = 6.0f;
        m_Param [ PGROUND_SLOPE ] = 0.10f;
        break;
    case 3:		// Small dam break
        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values( -40, 0, -40  );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values( 40, 60, 40 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( 0, 8, -35 );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 35, 55, 35 );
        m_Param [ PFORCE_MIN ] = 0.0f;
        m_Param [ PFORCE_MAX ] = 0.0f;
        m_Vec [ PPLANE_GRAV_DIR ] = Vector3DF_init_with_values( 0.0f, -9.8f, 0.0f );
        break;
    case 4:		// Dual-Wave pool
        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values( -100, 0, -15 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values( 100, 100, 15 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( -80, 8, -10 );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 80, 90, 10 );
        m_Param [ PFORCE_MIN ] = 20.0;
        m_Param [ PFORCE_MAX ] = 20.0;
        m_Vec [ PPLANE_GRAV_DIR ] = Vector3DF_init_with_values( 0.0f, -9.8f, 0.0f );
        break;
    case 5:		// Microgravity
        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values( -80, 0, -80 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values( 80, 100, 80 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( -60, 40, -60 );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 60, 80, 60 );
        m_Vec [ PPLANE_GRAV_DIR ] = Vector3DF_init_with_values( 0, -1, 0 );
        m_Param [ PGROUND_SLOPE ] = 0.1f;
        break;
    case 6:     // Morphogenesis small demo
        m_Param [ PSIMSCALE ] = 1.0f;
        m_Param [ PRADIUS ] = 1.0f;
        m_Param [ PSMOOTHRADIUS ] = 1.0f;
        m_Param [ PVISC ] = 0.1f;

        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values( 0, 0, 0 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values( 10, 20, 50 ); //( 80, 50, 80 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( m_Vec [ PVOLMIN ].x,  m_Vec [ PVOLMIN ].y, m_Vec [ PVOLMIN ].z );// will be reset to m_Vec[PBOUNDMIN].
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 10, 20, 30 );

        m_Param [ PGRAV ] = 2.000000f;
        m_Vec [ PPLANE_GRAV_DIR ] = Vector3DF_init_with_values( 0, -1, 0 );
        m_Param [ PGROUND_SLOPE ] = 0.1f;
        break;
    case 7:     // From SpecificationFile.txt
        m_Time = launchParams.m_Time;
        m_DT = launchParams.m_DT;
        m_Param [ PGRIDSIZE ] = launchParams.gridsize;
        m_Param [ PSPACING ] = launchParams.spacing;
        m_Param [ PSIMSCALE ] = launchParams.simscale;
        m_Param [ PSMOOTHRADIUS ] = launchParams.smoothradius;
        m_Param [ PVISC ] = launchParams.visc;
        m_Param [ PSURFACE_TENSION ] = launchParams.surface_tension;
        m_Param [ PMASS ] = launchParams.mass;
        m_Param [ PRADIUS ] = launchParams.radius;
        /*m_Param [ PDIST ] = launchParams.dist;*/
        m_Param [ PINTSTIFF ] = launchParams.intstiff;
        m_Param [ PEXTSTIFF ] = launchParams.extstiff;
        m_Param [ PEXTDAMP ] = launchParams.extdamp;
        m_Param [ PACCEL_LIMIT ] = launchParams.accel_limit;
        m_Param [ PVEL_LIMIT ] = launchParams.vel_limit;
        m_Param [ PGRAV ] = launchParams.grav;
        m_Param [ PGROUND_SLOPE ] = launchParams.ground_slope;
        m_Param [ PFORCE_MIN ] = launchParams.force_min;
        m_Param [ PFORCE_MAX ] = launchParams.force_max;
        m_Param [ PFORCE_FREQ ] = launchParams.force_freq;

        m_Vec [ PVOLMIN ] = launchParams.volmin;
        m_Vec [ PVOLMAX ] = launchParams.volmax;
        m_Vec [ PINITMIN ] = launchParams.initmin;
        m_Vec [ PINITMAX ] = launchParams.initmax;

        m_Param [ PACTUATION_FACTOR ] = launchParams.actuation_factor;
        m_Param [ PACTUATION_PERIOD ] = launchParams.actuation_period;

        break;
    case 8:  // default demo for parameter sweeps
        launchParams.num_particles = 4000;
        launchParams.demoType = 0;
        launchParams.simSpace = 7;

        m_Time = 0.000000;
        m_DT = 0.003000;

        m_Param [ PGRIDSIZE ] = 1.0;
        m_Param [ PSPACING ] = 1.000000;

        m_Param [ PSIMSCALE ] = 1.0;
        m_Param [ PSMOOTHRADIUS ] = 1.0;

        m_Param [ PVISC ] = 0.500000;
        m_Param [ PSURFACE_TENSION ] = 1.000000;

        m_Param [ PMASS ] = 0.00205;
        m_Param [ PRADIUS ] = 1.0000;

        m_Param [ PINTSTIFF ] = 2.000000;
        m_Param [ PEXTSTIFF ] = 50000.000000;
        m_Param [ PEXTDAMP ] = 100.000000;

        m_Param [ PACCEL_LIMIT ] = 150.000000;
        m_Param [ PVEL_LIMIT ] = 3.000000;

        m_Param [ PGRAV ] = 10.000000;
        m_Param [ PGROUND_SLOPE ] = 0.100000;

        m_Param [ PFORCE_MIN ] = 0.000000;
        m_Param [ PFORCE_MAX ] = 0.000000;
        m_Param [ PFORCE_FREQ ] = 16.000000;

        launchParams.x_dim = 10.000000;
        launchParams.y_dim = 10.000000;
        launchParams.z_dim = 3.000000;

        launchParams.pos_x = 0.000000;
        launchParams.pos_y = 0.000000;
        launchParams.pos_z = 0.000000;

        m_Vec [ PVOLMIN ] = Vector3DF_init_with_values( 0.0, 0.0, 0.0 );
        m_Vec [ PVOLMAX ] = Vector3DF_init_with_values( 10.0, 20.0, 20.0 );
        m_Vec [ PINITMIN ] = Vector3DF_init_with_values( 2.0, 2.0, 2.0 );
        m_Vec [ PINITMAX ] = Vector3DF_init_with_values( 10.0, 20.0, 10.0 );

        launchParams.num_files = 4000;
        launchParams.steps_per_InnerPhysicalLoop = 3;
        launchParams.steps_per_file = 6;
        launchParams.freeze_steps = 1;

        launchParams.debug = 0;
        launchParams.file_num = 0;

        launchParams.save_ply = 'n';
        launchParams.save_csv = 'n';
        launchParams.save_vtp = 'n';

        launchParams.gene_activity = 'n';
        launchParams.remodelling = 'n';
        launchParams.read_genome = 'y';

        m_Param [ PACTUATION_FACTOR ] = 0;
        m_Param [ PACTUATION_PERIOD ] = 1;

        break;
    }
    //std::cout<<"\nSetupExampleParams()3: launchParams.genomePath = "<<launchParams.genomePath<<"\n"<<std::flush;
}

void FluidSystem::SetupExampleGenome()  {   // need to set up a demo genome
    // Null genome
    for(int i=0; i< NUM_GENES; i++) m_FGenome.mutability[i] = 0;
    for(int i=0; i< NUM_GENES; i++) m_FGenome.delay[i] = 1;
    for(int i=0; i< NUM_GENES; i++) for(int j=0; j< NUM_GENES; j++) m_FGenome.sensitivity[i][j] = j;

    for(int i=0; i< NUM_TF/2; i++)      m_FGenome.tf_diffusability[i]    = 0;           // 1st half of TFs are non-diffusible.
    for(int i=NUM_TF/2; i< NUM_TF; i++) m_FGenome.tf_diffusability[i]    = 1;
    for(int i=0; i< NUM_TF; i++)        m_FGenome.tf_breakdown_rate[i]  = 1;

    for(int i=0; i< NUM_GENES; i++) for(int j=0; j< 2*NUM_TF+1; j++) m_FGenome.secrete[i][j]=0;     // 1st zero arrays.
    for(int i=0; i< NUM_GENES; i++) for(int j=0; j< 2*NUM_TF+1; j++) m_FGenome.activate[i][j]=0;

    m_FGenome.secrete[0][2*NUM_TF] = 2; // gene [0] secretes TF 1 & 3, at rates 1 & 4. // minimal test case.
    m_FGenome.secrete[0][2*0] = 1;
    m_FGenome.secrete[0][2*0+1] = 1;
    m_FGenome.secrete[0][2*1] = 3;
    m_FGenome.secrete[0][2*1+1] = 4;

    m_FGenome.activate[0][2*NUM_TF] = 1; // gene [0] activates TF 5, at rates 6. // minimal test case.
    m_FGenome.activate[0][2*0] = 5;
    m_FGenome.activate[0][2*0+1] = 6;

    //FBondParams *params_ =  &m_FGenome.fbondparams[0];
    //Particle remodelling, bond defaults & limits  m_FGenome.param[3][12]
    //0=elastin
    m_FGenome.param[0][m_FGenome.elongation_threshold]   = 0.5  ;
    m_FGenome.param[0][m_FGenome.elongation_factor]      = 0.002  ;
    m_FGenome.param[0][m_FGenome.strength_threshold]     = 1.0  ;
    m_FGenome.param[0][m_FGenome.strengthening_factor]   = 0.002  ;

    m_FGenome.param[0][m_FGenome.max_rest_length]        = 1.0  ;
    m_FGenome.param[0][m_FGenome.min_rest_length]        = 0.3  ;
    m_FGenome.param[0][m_FGenome.max_modulus]            = 100000;
    m_FGenome.param[0][m_FGenome.min_modulus]            = 10  ;

    m_FGenome.param[0][m_FGenome.elastLim]               = 8  ;
    m_FGenome.param[0][m_FGenome.default_rest_length]    = 0.5  ;
    m_FGenome.param[0][m_FGenome.default_modulus]        = 1000;
    m_FGenome.param[0][m_FGenome.default_damping]        = 10;

    //1=collagen
    m_FGenome.param[1][m_FGenome.elongation_threshold]   = 4.0  ;
    m_FGenome.param[1][m_FGenome.elongation_factor]      = 0.01 ;
    m_FGenome.param[1][m_FGenome.strength_threshold]     = 4.1  ;
    m_FGenome.param[1][m_FGenome.strengthening_factor]   = 0.01 ;

    m_FGenome.param[1][m_FGenome.max_rest_length]        = 1.0  ;
    m_FGenome.param[1][m_FGenome.min_rest_length]        = 0.3  ;
    m_FGenome.param[1][m_FGenome.max_modulus]            = 1000000000;//0.8  ;
    m_FGenome.param[1][m_FGenome.min_modulus]            = 1000000;//0.3  ;

    m_FGenome.param[1][m_FGenome.elastLim]               = 0.55  ;
    m_FGenome.param[1][m_FGenome.default_rest_length]    = 0.5  ;
    m_FGenome.param[1][m_FGenome.default_modulus]        = 10000000  ;
    m_FGenome.param[1][m_FGenome.default_damping]        = 100  ;

    //2=apatite
    m_FGenome.param[2][m_FGenome.elongation_threshold]   = 1.0  ;
    m_FGenome.param[2][m_FGenome.elongation_factor]      = 0.001  ;
    m_FGenome.param[2][m_FGenome.strength_threshold]     = 1.0  ;
    m_FGenome.param[2][m_FGenome.strengthening_factor]   = 0.001  ;

    m_FGenome.param[2][m_FGenome.max_rest_length]        = 1.0  ;
    m_FGenome.param[2][m_FGenome.min_rest_length]        = 0.3  ;
    m_FGenome.param[2][m_FGenome.max_modulus]            = 1000000000;//0.8  ;
    m_FGenome.param[2][m_FGenome.min_modulus]            = 1000000;//0.3  ;

    m_FGenome.param[2][m_FGenome.elastLim]               = 0.05  ;
    m_FGenome.param[2][m_FGenome.default_rest_length]    = 0.5  ;
    m_FGenome.param[2][m_FGenome.default_modulus]        = 10000000  ;
    m_FGenome.param[2][m_FGenome.default_damping]        = 1000  ;

    //Bond remodelling  m_FGenome.tanh_param[3][8];                 // lengthening/shortening
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.l_a] = 0.0;   // y-shift
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.l_b] = 0.0;   // y-scaling
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.l_c] = 0.0;   // x-scaling
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.l_d] = 0.0;   // x-shift

    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.s_a] = 1.008; // strengthening/weakening  // mod_mul = 1.008+0.01*np.tanh(10*si2[i+1]-7)
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.s_b] = 0.01;
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.s_c] = 10.0;
    m_FGenome.tanh_param[m_FGenome.elastin][m_FGenome.s_d] = 7.0;

    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.l_a] = 0.0;  // lengthening/shortening
    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.l_b] = 0.0;
    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.l_c] = 0.0;
    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.l_d] = 0.0;

    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.s_a] = 0.0;  // strengthening/weakening
    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.s_b] = 0.0;
    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.s_c] = 0.0;
    m_FGenome.tanh_param[m_FGenome.collagen][m_FGenome.s_d] = 0.0;

    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.l_a] = 0.0;   // lengthening/shortening
    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.l_b] = 0.0;
    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.l_c] = 0.0;
    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.l_d] = 0.0;

    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.s_a] = 0.0;   // strengthening/weakening
    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.s_b] = 0.0;
    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.s_c] = 0.0;
    m_FGenome.tanh_param[m_FGenome.apatite][m_FGenome.s_d] = 0.0;

}

void FluidSystem::SetupSpacing (){
    m_Param [ PSIMSIZE ] = m_Param [ PSIMSCALE ] * (m_Vec[PVOLMAX].z - m_Vec[PVOLMIN].z);

    if ( m_Param[PSPACING] == 0 ) {
        // Determine spacing from density
        m_Param [PDIST] = pow ( (float) m_Param[PMASS] / m_Param[PRESTDENSITY], 1/3.0f );
        m_Param [PSPACING] = m_Param [ PDIST ]*0.87f / m_Param[ PSIMSCALE ];
    } else {
        // Determine density from spacing
        m_Param [PDIST] = m_Param[PSPACING] * m_Param[PSIMSCALE] / 0.87f;
        m_Param [PRESTDENSITY] = m_Param[PMASS] / pow ( (float) m_Param[PDIST], 3.0f );
    }
    if (m_FParams.debug>0)printf ( "\nSetupSpacing: Density=,%f, Spacing=,%f, PDist=,%f\n", m_Param[PRESTDENSITY], m_Param[PSPACING], m_Param[PDIST] );

    // Particle Boundaries
    m_Vec[PBOUNDMIN] = m_Vec[PVOLMIN];
    m_Vec[PBOUNDMIN] = *Vector3DF_addDouble(&m_Vec[PBOUNDMIN],  2.0*(m_Param[PGRIDSIZE] / m_Param[PSIMSCALE])  );
    m_Vec[PBOUNDMAX] = m_Vec[PVOLMAX];
    m_Vec[PBOUNDMAX] = Vector3DF_subtractDouble(&m_Vec[PBOUNDMAX], 2.0*(m_Param[PGRIDSIZE] / m_Param[PSIMSCALE])  );
}

void FluidSystem::SetupSimulation(int gpu_mode, int cpu_mode){ // const char * relativePath, int gpu_mode, int cpu_mode             //FUNCTIONS COMMMENTED OUT IN HERE!!!!!!!!!!!!!!!!!
     // Allocate buffers for points
    //std::cout<<"\nSetupSimulation chk1, m_FParams.debug="<<m_FParams.debug<<std::flush;
    m_Param [PNUM] = launchParams.num_particles;                             // NB there is a line of text above the particles, hence -1.
    mMaxPoints = m_Param [PNUM];
    m_Param [PGRIDSIZE] = 2*m_Param[PSMOOTHRADIUS] / m_Param[PGRID_DENSITY];
    //std::cout<<"\nSetupSimulation chk2, m_FParams.debug="<<m_FParams.debug<<std::flush;

    SetupSPH_Kernels ();
    SetupSpacing ();
    SetupGrid ( m_Vec[PVOLMIN]/*bottom corner*/, m_Vec[PVOLMAX]/*top corner*/, m_Param[PSIMSCALE], m_Param[PGRIDSIZE]);
    //std::cout<<"\nSetupSimulation chk3, m_FParams.debug="<<m_FParams.debug<<std::flush;

    if (gpu_mode != GPU_OFF) {     // create CUDA instance etc..
        //FluidSetupCL ( mMaxPoints, m_GridSrch, *(int3*)& m_GridRes, *(float3*)& m_GridSize, *(float3*)& m_GridDelta, *(float3*)& m_GridMin, *(float3*)& m_GridMax, m_GridTotal, 0 );  CUCLCUCL
        UpdateParams();            //  sends simulation params to device.
        UpdateGenome();            //  sends genome to device.              // NB need to initialize genome from file, or something.
    }
    if (m_FParams.debug>1)std::cout<<"\nSetupSimulation chk4, mMaxPoints="<<mMaxPoints<<", gpu_mode="<<gpu_mode<<", cpu_mode="<<cpu_mode<<", m_FParams.debug="<<m_FParams.debug<<std::flush;

    AllocateParticles ( mMaxPoints, gpu_mode, cpu_mode );  // allocates only cpu buffer for particles
    if (m_FParams.debug>1)std::cout<<"\nSetupSimulation chk5 "<<std::flush;

    //AllocateGrid(gpu_mode, cpu_mode);
    //if (m_FParams.debug>1)std::cout<<"\nSetupSimulation chk6 "<<std::flush;

}
