#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 300
///////////////////////////////////////
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
//#include "RunCL.h"
#include "fluid.h"
#include "fluid_system.h"
//#include "fluid_system.cpp"
///////////////////////////////////////
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }
#define SDK_SUCCESS 0
#define SDK_FAILURE 1
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;


RunCL::RunCL(Json::Value obj_)
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

RunCL::~RunCL()
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

void RunCL::allocatemem(FParams fparam, FBufs *fbuf, FBufs *ftemp, FGenome fgenome)
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


																				if(verbosity>0) cout << "RunCL::allocatemem_finished\n\n" << flush;
}


void RunCL::CleanUp()
{
																																			cout<<"\nRunCL::CleanUp_chk0"<<flush;
	cl_int status;
	status = clReleaseMemObject(m_FParamDevice);	if (status != CL_SUCCESS)	{ cout << "\nbasemem  status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.1"<<flush;
	status = clReleaseMemObject(m_FluidDevice);	if (status != CL_SUCCESS)	{ cout << "\nimgmem   status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.2"<<flush;
	status = clReleaseMemObject(m_FluidTempDevice);	if (status != CL_SUCCESS)	{ cout << "\ncdatabuf status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.3"<<flush;
	status = clReleaseMemObject(m_FGenomeDevice);	if (status != CL_SUCCESS)	{ cout << "\nhdatabuf status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.4"<<flush;
																																			cout<<"\nRunCL::CleanUp_chk1_finished"<<flush;
}

void RunCL::exit_(cl_int res)
{
	CleanUp();
	//~RunCL(); Never need to call a destructor manually.
	exit(res);
}

