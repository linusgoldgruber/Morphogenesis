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
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
#include "fluid.h"
#include "fluid_system.h"
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }
//#include <curand_kernel.h> //NOT REPLACED YET. POSSIBLE OPTIONS: https://acesse.dev/Y9Zcq

#define SDK_SUCCESS 0
#define SDK_FAILURE 1
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

void FluidSystem::allocatemem(FParams fparam, FBufs *fbuf, FBufs *ftemp, FGenome fgenome)
{
																				if(verbosity>0) cout << "RunCL::allocatemem_chk0" << flush;
	stringstream 		ss;
	ss << "allocatemem";
	cl_int status;
	cl_event writeEvt;

	int layerstep 		= fparam.szPnts * height;
	// Get the maximum work group size for executing the kernel on the device ///////// From https://github.com/rsnemmen/OpenCL-examples/blob/e2c34f1dfefbd265cfb607c2dd6c82c799eb322a/square_array/square.c
	status = clGetKernelWorkGroupInfo(m_Kern[FUNC_FPREFIXUP], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clGetKernelWorkGroupInfo(m_Kern[FUNC_FPREFIXSUM], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clGetKernelWorkGroupInfo(m_Kern[FUNC_TALLYLISTS], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}

	// Number of total work items, calculated here after 1st image is loaded &=> know the size.
	// NB localSize must be devisor
	// NB global_work_size must be a whole number of "Preferred work group size multiple" for Nvidia.
	global_work_size = ceil((float)layerstep/(float)local_work_size) * local_work_size;
	local_work_size=32; // trial for nvidia
																				if(verbosity>1){
																					cout<<"\nglobal_work_size="<<global_work_size<<", local_work_size="<<local_work_size<<", deviceId="<<deviceId<<"\n"<<flush;
																					cout<<"\nlayerstep=width*height="<<width<<"*"<<height<<"="<<layerstep<<",\tsizeof(layerstep)="<< sizeof(layerstep) <<",\tsizeof(int)="<< sizeof(int) <<flush;
																					cout<<"\n";
																				}
	cl_int res;
	m_FParamDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(fparam), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	m_FluidDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(fbuf), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	m_FluidTempDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(ftemp), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	m_FGenomeDevice		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , sizeof(fgenome), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}

																				if(verbosity>1) {
																					cout << ",m_FParamDevice = " 		<< m_FParamDevice << endl;
																					cout << ",m_FluidDevice = " 		<< m_FluidDevice << endl;
																					cout << ",m_FluidTempDevice = " 	<< m_FluidTempDevice << endl;
																					cout << ",m_FGenomeDevice = " 	<< m_FGenomeDevice << endl;
																				}

	status = clEnqueueWriteBuffer(uload_queue, m_FParamDevice, 		CL_FALSE, 0, sizeof(fparam), 		&fparam, 			0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.3\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, m_FluidDevice, 		CL_FALSE, 0, sizeof(fbuf), 		&fbuf, 			0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.4\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, m_FluidTempDevice, 	CL_FALSE, 0, sizeof(ftemp), 				&ftemp, 		0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.5\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, m_FGenomeDevice, 	CL_FALSE, 0, sizeof(fgenome)*3, 		&fgenome, 0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.6\n" << endl;exit_(status);}

	clFlush(uload_queue); status = clFinish(uload_queue); 					if (status != CL_SUCCESS)	{ cout << "\nclFinish(uload_queue)="				<< status << checkerror(status)<<"\n"<<flush; exit_(status);}

	// set kernelArg. NB "0 &k2kbuf" & "2 &imgmem" set in calcCostVol(..)
	res = clSetKernelArg(m_Kern[FUNC_INSERT], 1, sizeof(cl_mem),  &m_FParamDevice);		if(res!=CL_SUCCESS){cout<<"\nfparam res= "   		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	cout << "TEST: "<< checkerror(res);
	res = clSetKernelArg(m_Kern[FUNC_FPREFIXUP], 2, sizeof(cl_mem),  &m_FluidDevice);		if(res!=CL_SUCCESS){cout<<"\nfbuf res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	res = clSetKernelArg(m_Kern[FUNC_TALLYLISTS], 1, sizeof(cl_mem),  &m_FParamDevice);		if(res!=CL_SUCCESS){cout<<"\nfparam res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	res = clSetKernelArg(m_Kern[FUNC_TALLYLISTS], 2, sizeof(cl_mem),  &m_FluidDevice);		if(res!=CL_SUCCESS){cout<<"\nfbuf res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	res = clSetKernelArg(m_Kern[FUNC_COUNTING_SORT], 2, sizeof(cl_mem),  &m_FParamDevice);		if(res!=CL_SUCCESS){cout<<"\nfparam res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	res = clSetKernelArg(m_Kern[FUNC_COUNTING_SORT], 2, sizeof(cl_mem),  &m_FluidDevice);		if(res!=CL_SUCCESS){cout<<"\nfbuf res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	res = clSetKernelArg(m_Kern[FUNC_COUNTING_SORT], 2, sizeof(cl_mem),  &m_FluidTempDevice);		if(res!=CL_SUCCESS){cout<<"\nftemp res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //
	res = clSetKernelArg(m_Kern[FUNC_CONTRIBUTE_PRESSURE], 4, sizeof(cl_mem),  &m_FluidTempDevice);		if(res!=CL_SUCCESS){cout<<"\nftemp res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} //


																				if(verbosity>0) cout << "RunCL::allocatemem_finished\n\n" << flush;
}

void FluidSystem::CleanUp()
{
																																			cout<<"\nRunCL::CleanUp_chk0"<<flush;
	cl_int status;
	status = clReleaseMemObject(m_FParamDevice);	if (status != CL_SUCCESS)	{ cout << "\nbasemem  status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.1"<<flush;
	status = clReleaseMemObject(m_FluidDevice);	if (status != CL_SUCCESS)	{ cout << "\nimgmem   status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.2"<<flush;
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
            clEnqueueWriteBuffer(queue, dest_buf, CL_TRUE, 0, cpucnt*stride, src_buf, 0, NULL, NULL);
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
