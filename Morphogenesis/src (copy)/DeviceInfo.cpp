/*
 *
 *
 *
*/

#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 300

#include <CL/opencl.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "/home/goldi/Documents/KDevelop Projects/HandsOnOpenCL/Exercises-Solutions/Exercises/Cpp_common/err_code.h"
#include <stdio.h>
#include <regex>
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }

int main() {
    cl_int err;
    cl_uint num_platforms;
    cl_platform_id platform;
    cl_device_id gpu_device;
    cl_context context;
    cl_command_queue gpu_queue;
    cl_program program;
    cl_kernel kernel;
    cl_mem input_buffer;
    int input_data[16] = {0};


    // Read the kernel source from the file
    std::ifstream kernelFile("/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/fluid_system_opencl.cl");
    std::string kernelSource((std::istreambuf_iterator<char>(kernelFile)),
                             std::istreambuf_iterator<char>());

    // Output the size of the kernel source file
    std::cout << "Size of kernel source file: " << kernelSource.size() << " bytes" << std::endl;

    // Count the number of kernels defined in the source
    std::regex kernelRegex("\\b__kernel\\b");
    auto kernelBegin = std::sregex_iterator(kernelSource.begin(), kernelSource.end(), kernelRegex);
    auto kernelEnd = std::sregex_iterator();
    int kernelCount = std::distance(kernelBegin, kernelEnd);

    // Output the number of kernels
    std::cout << "Number of kernels: " << kernelCount << std::endl;

    const char* kernelSourceStr = kernelSource.c_str();

   // Get platform
    err = clGetPlatformIDs(1, &platform, &num_platforms);
    CHECK_ERROR(err);

    // Get GPU device
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &gpu_device, NULL);
    CHECK_ERROR(err);

    // Create context
    context = clCreateContext(NULL, 1, &gpu_device, NULL, NULL, &err);
    CHECK_ERROR(err);

    // Create command queue
    cl_queue_properties properties[] = {0};
    gpu_queue = clCreateCommandQueueWithProperties(context, gpu_device, properties, &err);
    CHECK_ERROR(err);



    // Create program
    program = clCreateProgramWithSource(context, 1, &kernelSourceStr, NULL, &err);
    CHECK_ERROR(err);

        // Output the number of kernels
    //std::cout << "Kernels: \n" << kernelSourceStr << std::endl;

    // Build program
    err = clBuildProgram(program, 1, &gpu_device, "", NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t log_size;
        char *log;
        clGetProgramBuildInfo(program, gpu_device, CL_PROGRAM_BUILD_LOG,
                              0, NULL,&log_size);
        log = (char*)malloc(log_size+1);
        clGetProgramBuildInfo(program, gpu_device,
                              CL_PROGRAM_BUILD_LOG,
                              log_size+1,
                              log,
                              NULL);
        printf("%s\n", log);
        free(log);
        exit(1);
    }

    // Get maximum work-group size for device
    size_t max_work_group_size;
    clGetDeviceInfo(gpu_device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, NULL);
    double work_size = max_work_group_size;
    printf("Max size is %f",work_size);

    // Create kernel
    kernel = clCreateKernel(program,"increment",&err);
    CHECK_ERROR(err);

   // Debug Output
    std::cout << "After Kernel\n " << std::endl;

    // Create buffer
    input_buffer = clCreateBuffer(context,
                                  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                  sizeof(input_data),
                                  input_data,
                                  &err);
    CHECK_ERROR(err);

   // Set kernel arguments
   err = clSetKernelArg(kernel,
                        0,
                        sizeof(cl_mem),
                        &input_buffer);
   CHECK_ERROR(err);

   // Enqueue kernel to queue
   size_t global_size = 16;
   err = clEnqueueNDRangeKernel(gpu_queue,
                                kernel,
                                1,
                                NULL,
                                &global_size,
                                NULL,
                                0,
                                NULL,
                                NULL);
   CHECK_ERROR(err);

   // Wait for queue to finish
   err = clFinish(gpu_queue);
   CHECK_ERROR(err);

   // Read results
   err = clEnqueueReadBuffer(gpu_queue,
                             input_buffer,
                             CL_TRUE,
                             0,
                             sizeof(input_data),
                             input_data,
                             0,
                             NULL,
                             NULL);
   CHECK_ERROR(err);

   for (int i=0; i<16; i++) {
       printf("%d ", input_data[i]);
   }
   printf("\n");

   // Debug Output
    std::cout << "End of File\n " << std::endl;

   return 0;
}
