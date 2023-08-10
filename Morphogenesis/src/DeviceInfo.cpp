#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 300

#include <CL/opencl.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <regex>
#include <jsoncpp/json/json.h>
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }

#include "RunCL.h"

#define SDK_SUCCESS 0
#define SDK_FAILURE 1

using namespace std;

RunCL::RunCL(Json::Value obj_)
{
	obj = obj_;
	verbosity = obj["verbosity"].asInt();
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
	m_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);		if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);}
	uload_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);	if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);}
	dload_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);	if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);}
	track_queue = clCreateCommandQueueWithProperties(m_context, deviceId, prop, &status);	if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);}
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

	const char includeOptions[] = "-I /lib/i386-linux-gnu"; // Paths to include directories ///// -I /usr/lib/gcc/x86_64-linux-gnu/11/include /////

	status = clBuildProgram(m_program, 1, devices, includeOptions, NULL, NULL);					/*Step 6: Build program.*/////////////////////
	if (status != CL_SUCCESS){
		printf("\nclBuildProgram failed: %d\n", status);
		char buf[0x10000];
		clGetProgramBuildInfo(m_program, deviceId, CL_PROGRAM_BUILD_LOG, 0x10000, buf, NULL);
		printf("\n%s\n End of clBuildProgram error log..", buf);
		exit_(status);
	}

	insertParticlesCL_kernel     = clCreateKernel(m_program, "insertParticlesCLTest", NULL);				/*Step 7: Create kernel objects.*////////////
	prefixFixup_kernel     = clCreateKernel(m_program, "prefixFixupTest", NULL);

	basemem=imgmem=cdatabuf=hdatabuf=k2kbuf=dmem=amem=basegraymem=gxmem=gymem=g1mem=lomem=himem=0;		// set device pointers to zero
																						if(verbosity>0) cout << "RunCL_constructor finished\n" << flush;
}


RunCL::~RunCL()
{
	cl_int status;
	status = clReleaseKernel(prefixFixup_kernel);      	if (status != CL_SUCCESS)	{ cout << "\nRelease Kernel1 status = " << checkerror(status) <<"\n"<<flush; }
	//status = clReleaseKernel(cache3_kernel);		if (status != CL_SUCCESS)	{ cout << "\nRelease Kernel2 status = " << checkerror(status) <<"\n"<<flush; }

	status = clReleaseProgram(m_program);			if (status != CL_SUCCESS)	{ cout << "\nRelease Program status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(m_queue);		if (status != CL_SUCCESS)	{ cout << "\nRelease CQ1 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(uload_queue);	if (status != CL_SUCCESS)	{ cout << "\nRelease CQ2 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(dload_queue);	if (status != CL_SUCCESS)	{ cout << "\nRelease CQ3 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(track_queue);	if (status != CL_SUCCESS)	{ cout << "\nRelease CQ4 status = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseContext(m_context);			if (status != CL_SUCCESS)	{ cout << "\nRelease Context status = " << checkerror(status) <<"\n"<<flush; }
}

void RunCL::CleanUp()
{
																																			cout<<"\nRunCL::CleanUp_chk0"<<flush;
	cl_int status;
	status = clReleaseMemObject(basemem);	if (status != CL_SUCCESS)	{ cout << "\nbasemem  status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.1"<<flush;
	status = clReleaseMemObject(imgmem);	if (status != CL_SUCCESS)	{ cout << "\nimgmem   status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.2"<<flush;
	status = clReleaseMemObject(cdatabuf);	if (status != CL_SUCCESS)	{ cout << "\ncdatabuf status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.3"<<flush;
	status = clReleaseMemObject(hdatabuf);	if (status != CL_SUCCESS)	{ cout << "\nhdatabuf status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.4"<<flush;
	status = clReleaseMemObject(k2kbuf);	if (status != CL_SUCCESS)	{ cout << "\nk2kbuf   status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.5"<<flush;
	status = clReleaseMemObject(qmem);		if (status != CL_SUCCESS)	{ cout << "\ndmem     status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.7"<<flush;
	status = clReleaseMemObject(dmem);		if (status != CL_SUCCESS)	{ cout << "\ndmem     status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.8"<<flush;
	status = clReleaseMemObject(amem);		if (status != CL_SUCCESS)	{ cout << "\namem     status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.9"<<flush;
	status = clReleaseMemObject(lomem);		if (status != CL_SUCCESS)	{ cout << "\nlomem    status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.10"<<flush;
	status = clReleaseMemObject(himem);		if (status != CL_SUCCESS)	{ cout << "\nhimem    status = " << checkerror(status) <<"\n"<<flush; }		if(verbosity>0) cout<<"\nRunCL::CleanUp_chk0.11"<<flush;
																																			cout<<"\nRunCL::CleanUp_chk1_finished"<<flush;
}

void RunCL::exit_(cl_int res)
{
	CleanUp();
	//~RunCL(); Never need to call a destructor manually.
	exit(res);
}


int main()
{
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;
    obj["kernel_filepath"] = "/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/kernelTest.cl";
    RunCL runCLInstance(obj);

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// int main() {
//     cl_uint numPlatforms;
//     clGetPlatformIDs(0, NULL, &numPlatforms);
//     std::cout << "Number of platforms: " << numPlatforms << std::endl;
//
//     cl_platform_id* platforms = new cl_platform_id[numPlatforms];
//     clGetPlatformIDs(numPlatforms, platforms, NULL);
//
//     for (cl_uint i = 0; i < numPlatforms; i++) {
//         char platformName[128];
//         clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(platformName), platformName, NULL);
//         std::cout << "Platform " << i << ": " << platformName << std::endl;
//
//         cl_uint numDevices;
//         clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
//         std::cout << "  Number of devices: " << numDevices << std::endl;
//
//         cl_device_id* devices = new cl_device_id[numDevices];
//         clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
//
//         for (cl_uint j = 0; j < numDevices; j++) {
//             char deviceName[128];
//             clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
//             std::cout << "  Device " << j << ": " << deviceName << std::endl;
//         }
//
//         delete[] devices;
//     }
//
//     delete[] platforms;
//     return 0;
// }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//     cl_int err;
//     cl_uint num_platforms;
//     cl_platform_id platform;
//     cl_device_id gpu_device;
//     cl_context context;
//     cl_command_queue gpu_queue;
//     cl_program program;
//     cl_kernel kernel;
//     cl_mem input_buffer;
//     int input_data[16] = {0};
//
//
//     // Read the kernel source from the file
//     std::ifstream kernelFile("/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/fluid_system_opencl.cl");
//     std::string kernelSource((std::istreambuf_iterator<char>(kernelFile)),
//                              std::istreambuf_iterator<char>());
//
//     // Output the size of the kernel source file
//     std::cout << "Size of kernel source file: " << kernelSource.size() << " bytes" << std::endl;
//
//     // Count the number of kernels defined in the source
//     std::regex kernelRegex("\\b__kernel\\b");
//     auto kernelBegin = std::sregex_iterator(kernelSource.begin(), kernelSource.end(), kernelRegex);
//     auto kernelEnd = std::sregex_iterator();
//     int kernelCount = std::distance(kernelBegin, kernelEnd);
//
//     // Output the number of kernels
//     std::cout << "Number of kernels: " << kernelCount << std::endl;
//
//     const char* kernelSourceStr = kernelSource.c_str();
//
//    // Get platform
//     err = clGetPlatformIDs(1, &platform, &num_platforms);
//     CHECK_ERROR(err);
//
//     // Get GPU device
//     err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &gpu_device, NULL);
//     CHECK_ERROR(err);
//
//     // Create context
//     context = clCreateContext(NULL, 1, &gpu_device, NULL, NULL, &err);
//     CHECK_ERROR(err);
//
//     // Create command queue
//     cl_queue_properties properties[] = {0};
//     gpu_queue = clCreateCommandQueueWithProperties(context, gpu_device, properties, &err);
//     CHECK_ERROR(err);
//
//
//
//     // Create program
//     program = clCreateProgramWithSource(context, 1, &kernelSourceStr, NULL, &err);
//     CHECK_ERROR(err);
//
//         // Output the number of kernels
//     //std::cout << "Kernels: \n" << kernelSourceStr << std::endl;
//
//     // Build program
//     err = clBuildProgram(program, 1, &gpu_device, "", NULL, NULL);
//     if (err != CL_SUCCESS) {
//         size_t log_size;
//         char *log;
//         clGetProgramBuildInfo(program, gpu_device, CL_PROGRAM_BUILD_LOG,
//                               0, NULL,&log_size);
//         log = (char*)malloc(log_size+1);
//         clGetProgramBuildInfo(program, gpu_device,
//                               CL_PROGRAM_BUILD_LOG,
//                               log_size+1,
//                               log,
//                               NULL);
//         printf("%s\n", log);
//         free(log);
//         exit(1);
//     }
//
//     // Get maximum work-group size for device
//     size_t max_work_group_size;
//     clGetDeviceInfo(gpu_device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, NULL);
//     double work_size = max_work_group_size;
//     printf("Max size is %f",work_size);
//
//     // Create kernel
//     kernel = clCreateKernel(program,"increment",&err);
//     CHECK_ERROR(err);
//
//    // Debug Output
//     std::cout << "After Kernel\n " << std::endl;
//
//     // Create buffer
//     input_buffer = clCreateBuffer(context,
//                                   CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                                   sizeof(input_data),
//                                   input_data,
//                                   &err);
//     CHECK_ERROR(err);
//
//    // Set kernel arguments
//    err = clSetKernelArg(kernel,
//                         0,
//                         sizeof(cl_mem),
//                         &input_buffer);
//    CHECK_ERROR(err);
//
//    // Enqueue kernel to queue
//    size_t global_size = 16;
//    err = clEnqueueNDRangeKernel(gpu_queue,
//                                 kernel,
//                                 1,
//                                 NULL,
//                                 &global_size,
//                                 NULL,
//                                 0,
//                                 NULL,
//                                 NULL);
//    CHECK_ERROR(err);
//
//    // Wait for queue to finish
//    err = clFinish(gpu_queue);
//    CHECK_ERROR(err);
//
//    // Read results
//    err = clEnqueueReadBuffer(gpu_queue,
//                              input_buffer,
//                              CL_TRUE,
//                              0,
//                              sizeof(input_data),
//                              input_data,
//                              0,
//                              NULL,
//                              NULL);
//    CHECK_ERROR(err);
//
//    for (int i=0; i<16; i++) {
//        printf("%d ", input_data[i]);
//    }
//    printf("\n");

   // Debug Output
//     std::cout << "End of File\n " << std::endl;
//
//    return 0;
