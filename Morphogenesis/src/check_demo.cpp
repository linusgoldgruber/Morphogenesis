#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 300
#define SDK_SUCCESS 0
#define SDK_FAILURE 1
#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <inttypes.h>
#include <errno.h>
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
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

    char paramsPath[256];
    char genomePath[256];
    char pointsPath[256];
    char outPath[256];

	if ( argc != 4 ){
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
    ifstream ifs(argv[3]);
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;
    bool b = reader.parse(ifs, obj);
    if (!b) { cout << "Error: " << reader.getFormattedErrorMessages();}   //else {cout << "NB lists .json file entries alphabetically: \n" << obj ;}
    cout << "\n\n\n" << endl;

    FluidSystem fluid(obj);

    uint debug = 2;
    fluid.SetDebug ( debug );

    // Clear all buffers
    fluid.Initialize();   // where do the buffers for params and genome get allocated when there is no fluid.InitializeOpenCL (); ?

    fluid.ReadSimParams(paramsPath);
    fluid.ReadGenome(genomePath);


    //fluid.ReadPointsCSV2_DEBUG(pointsPath, GPU_OFF, CPU_YES);  //fluid.ReadPointsCSV(pointsPath, GPU_OFF, CPU_YES);
    fluid.ReadPointsCSV2(pointsPath, GPU_SINGLE, CPU_YES);  //fluid.ReadPointsCSV(pointsPath, GPU_OFF, CPU_YES);
    printf("\nchk1: ReadPointsCSV2() \n");
    fluid.WriteSimParams ( outPath );
    printf("\nchk2: WriteSimParams() \n");
    fluid.WriteGenome( outPath );
    printf("\nchk3: WriteGenome() \n");
    fluid.SavePointsCSV2 ( outPath, 1 );
    printf("\nchk4: SavePointsCSV2() \n");
    fluid.SavePointsVTP2(outPath, 1 ); //TODO
    printf("\nchk5: SavePointsVTP2() \n");

    printf("\ncheck_demo finished.\n");
    fluid.Exit_no_CL ();
    return 0;

}
