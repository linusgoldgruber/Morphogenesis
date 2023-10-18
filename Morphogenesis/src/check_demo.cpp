#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 300
#define SDK_SUCCESS 0
#define SDK_FAILURE 1

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <vector_types.h>
#include <inttypes.h>
#include <errno.h>
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }
#include "fluid.h"
#include "fluid_system.h"
#include "host_CL.cpp"
#include <CL/cl.h>
#include <CL/opencl.h>
#include <chrono>
#include <filesystem>


//#include "fluid_system.h"

int main ( int argc, const char** argv ) 
{
    uint debug =2;
    char paramsPath[256];
    char genomePath[256];
    char pointsPath[256];
    char outPath[256];
	if ( argc != 3 ){
	    printf("usage: check_demo  simulation_data_folder  output_folder\n");
	    return 0;
	}else {
        sprintf ( paramsPath, "%s/SimParams.txt", argv[1] );
        printf("simulation parameters file = %s\n", paramsPath);
        
        sprintf ( genomePath, "%s/genome.csv", argv[1] );
        printf("simulation parameters file = %s\n", genomePath);
        
        sprintf ( pointsPath, "%s/particles_pos_vel_color100001.csv", argv[1] );  // particles_pos_vel_color100001_test_data.csv
        printf("simulation points file = %s\n", pointsPath);
        
        sprintf ( outPath, "%s", argv[2] );
        printf("output_folder = %s\n", outPath);
	}	
    
    // Initialize
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;
    obj["kernel_filepath"] = "/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/kernelTest.cl";
    FluidSystem fluid(obj);

    // Clear all buffers
    fluid.Initialize();   // where do the buffers for params and genome get allocated when there is no fluid.InitializeOpenCL (); ?

    fluid.ReadSimParams(paramsPath);

    fluid.ReadGenome(genomePath);

    cout << "##################################SO FAR So GOOD########################################";

    fluid.ReadPointsCSV2_TEST(pointsPath, GPU_OFF, CPU_YES);  //fluid.ReadPointsCSV(pointsPath, GPU_OFF, CPU_YES);
    printf("\nchk1\n");
    fluid.WriteSimParams ( outPath );
    printf("\nchk2\n");
    fluid.WriteGenome( outPath );
    printf("\nchk3\n");
    fluid.SavePointsCSV2 ( outPath, 1 );
    printf("\nchk4\n");
    fluid.SavePointsVTP2(outPath, 1 );
    printf("\nchk5\n");
    
    printf("\ncheck_demo finished.\n");
    //fluid.Exit_no_CL ();
    return 0;


    /*// Function
    void allocateAndAccessMemory(FBufs* fb, int n) {

        // Allocate with malloc
        fb->mcpu[n] = (char*)malloc(sizeof(char) * 20);

        // Check if memory allocation was successful
        if (fb->mcpu[n] != NULL) {

            // Write data
            strcpy(fb->mcpu[n], "Hello, dynamic memory!");

            // Read and print data
            std::cout << "Read from memroy: " << fb->mcpu[n] << std::endl;

            // Clean up (needed?)
            free(fb->mcpu[n]);
            fb->mcpu[n] = NULL;
        } else {
            std::cerr << "alloc failed." << std::endl;
        }
    }

int main ()
{




        FBufs fb;
        int n = 0; // Index of the memory buffer

        // Call the memory allocation and access function
        allocateAndAccessMemory(&fb, n);

        return 0;
        */
}
