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
#define CHECK_ERROR(status) if (status != CL_SUCCESS) { printf("Error: %d\n", status); exit(1); }
// #include "fluid.h"
#include "fluid_system.h"
#include "host_CL.cpp"
#include <chrono>
#include <filesystem>

int main ( int argc, const char** argv )
{
    char input_folder[256];
    char output_folder[256];
    char json_folder[256];
    cout << "argc: " << argc << "\n" << flush;

    if ((argc != 4) && (argc !=3)) {
        printf ( "usage: make_demo2 input_folder output_folder.\
        \nNB input_folder must contain \"SpecificationFile.txt\", output will be wrtitten to \"output_folder/out_data_time/\".\
        \nIf output_folder is not given the value from SpecificationFile.txt will be used.\n" );
        return 0;
    } else {
        sprintf ( input_folder, "%s", argv[1] );
        sprintf ( output_folder, "%s", argv[2] );
        sprintf ( json_folder, "%s", argv[3] );

        printf ( "input_folder = %s , output_folder = %s, \njson_folder = %s\n", input_folder, output_folder, json_folder );

    }

    // Initialize
    ifstream ifs(argv[3]);
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;

//     string line;
//     while (getline(ifs, line)) {
//         cout << line << endl;
//     }

    bool b = reader.parse(ifs, obj);
    if (!b) { cout << "Error: " << reader.getFormattedErrorMessages();}   else {cout << "NB lists .json file entries alphabetically: \n" << obj ;}
    cout << "\n\n\n" << endl;
    //obj["kernel_filepath"] = "/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/kernelTest.cl";
    FluidSystem fluid(obj);
    fluid.Initialize();
    fluid.InitializeOpenCL();
    uint debug =2;

    //std::cout<<"\n\nmake_demo2 chk0,"<<std::flush;

    fluid.ReadSpecificationFile ( input_folder );
    fluid.launchParams.debug = debug;
    std::cout<<"\n\nmake_demo2 chk1, fluid.launchParams.debug="<<fluid.launchParams.debug<<", fluid.launchParams.genomePath=" <<fluid.launchParams.genomePath  << ",  fluid.launchParams.spacing="<<fluid.launchParams.spacing<< ", fluid.launchParams.paramsPath="<<fluid.launchParams.paramsPath<< ",  fluid.launchParams.num_particles="<<fluid.launchParams.num_particles<< ",  fluid.launchParams.demoType="<<fluid.launchParams.demoType<<std::flush;

    for(int i=0; i<256; i++){fluid.launchParams.paramsPath[i] = input_folder[i];}
    for(int i=0; i<256; i++){fluid.launchParams.pointsPath[i] = input_folder[i];}
    //for(int i=0; i<256; i++){fluid.launchParams.genomePath[i] = input_folder[i];} // obtained from SpecificationFile.txt above.
    if(argc==4 || argc == 3){
        for(int i=0; i<256; i++) fluid.launchParams.outPath[i] = output_folder[i];
        cout << "\nfluid.launchParams.outPath: " << fluid.launchParams.outPath << flush;
    }
    if(mkdir(output_folder, 0755) == -1) cerr << "\nError :  failed to create output_folder.\n" << strerror(errno) << endl;
    else cout << "output_folder created\n"; // NB 0755 = rwx owner, rx for others.

    fluid.WriteDemoSimParams(           // Generates the simulation from data previously loaded from SpecificationFile.txt .
        fluid.launchParams.paramsPath, GPU_SINGLE, CPU_YES, fluid.launchParams.num_particles, fluid.launchParams.spacing, fluid.launchParams.x_dim, fluid.launchParams.y_dim, fluid.launchParams.z_dim, fluid.launchParams.demoType, fluid.launchParams.simSpace, fluid.launchParams.debug
    ); /*const char * relativePath*/
    //std::cout<<"\n\nmake_demo2 chk2 "<<std::flush;
    uint num_particles_start=fluid.ActivePoints();

    fluid.CreateFluidBuffers();
    fluid.TransferToCL ();
    fluid.Run2Simulation ();

    //std::cout<<"\n\nmake_demo2 chk3 "<<std::flush;
    fluid.WriteResultsCSV(input_folder, output_folder, num_particles_start);// NB post-slurm script to (i) cat results.csv files, (ii)tar-gzip and ftp folders to recipient.

    printf ( "\nClosed createDemo2.\n" );
    return 0;
}
