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

int main ( int argc, const char** argv ) 
{

    char paramsPath[256];
    char genomePath[256];
    char pointsPath[256];
    char outputFolder[256];

    //Find working directory-----------------------------
    fs::path Path = fs::current_path().parent_path();
    const char* directory = Path.c_str();
    std::cout << "Directory: " << directory << std::endl;
    //---------------------------------------------------

    //Initialize JSON---------------------------------------------------
    char jsonPath[512]; sprintf(jsonPath, "%s/%s", directory, argv[1]);
    cout << "Concatenated path: " << jsonPath << std::endl;
    ifstream ifs(jsonPath);
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;
    bool b = reader.parse(ifs, obj);
    if (!b) {cout << "Error: " << reader.getFormattedErrorMessages();}
    else    {cout << "NB lists .json file entries alphabetically: \n" << obj ;}
    cout << "\n\n\n" << endl;
    const char* in_path = obj["demo_path"].asCString();
    std::string out_path_str = obj.isMember("demo_path") && obj["demo_path"].isString() ? std::string(obj["demo_path"].asCString()) + "/check" : (std::cerr << "Error: 'demo_path' is missing or invalid." << std::endl, nullptr);
    const char* out_path = out_path_str.c_str();


    if ( argc != 2 ){
	    printf("usage: check_demo  JSON.config\n");
	    return 0;
	}else {
        sprintf ( paramsPath, "%s/%s/SimParams.txt", directory, in_path );
        printf("simulation parameters file = %s\n", paramsPath);

        sprintf ( genomePath, "%s/%s/genome.csv", directory, in_path );
        printf("simulation parameters file = %s\n", genomePath);

        sprintf ( pointsPath, "%s/%s/particles_pos_vel_color100001.csv", directory, in_path );  // particles_pos_vel_color100001_test_data.csv
        printf("simulation points file = %s\n", pointsPath);

        sprintf ( outputFolder, "%s/%s", directory, out_path );
        printf("outputFolder = %s\n", outputFolder);
	}

    // Directory setup
    if (!fs::exists(outputFolder)) {
            if (!fs::create_directories(outputFolder)) {std::cerr << "Error: Failed to create output folder: " << outputFolder << std::endl;return 1;}}



    FluidSystem fluid(obj);

    uint debug = 2;
    fluid.SetDebug ( debug );

    // Clear all buffers
    fluid.Initialize();   // where do the buffers for params and genome get allocated when there is no fluid.InitializeOpenCL (); ?
    fluid.InitializeOpenCL();

    fluid.ReadSimParams(paramsPath);
    fluid.ReadGenome(genomePath);


    //fluid.ReadPointsCSV2_DEBUG(pointsPath, GPU_OFF, CPU_YES);  //fluid.ReadPointsCSV(pointsPath, GPU_OFF, CPU_YES);
    fluid.ReadPointsCSV2(pointsPath, GPU_SINGLE, CPU_YES);  //fluid.ReadPointsCSV(pointsPath, GPU_OFF, CPU_YES);
    printf("\nchk1: ReadPointsCSV2() \n");
    fluid.WriteSimParams ( outputFolder );
    printf("\nchk2: WriteSimParams() \n");
    fluid.WriteGenome( outputFolder );
    printf("\nchk3: WriteGenome() \n");
    fluid.SavePointsCSV2 ( outputFolder, 1 );
    printf("\nchk4: SavePointsCSV2() \n");
    fluid.SavePointsVTP2(outputFolder, 1 ); //TODO
    printf("\nchk5: SavePointsVTP2() \n");

    printf("\ncheck_demo finished.\n");
    fluid.Exit_no_CL ();
    return 0;

}
