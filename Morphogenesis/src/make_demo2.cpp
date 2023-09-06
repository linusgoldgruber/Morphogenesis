#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include <string.h>
#include <chrono>
#include <filesystem>

#include "fluid_system.h"

typedef	unsigned int		uint;	

int main ( int argc, const char** argv ) 
{
    char input_folder[256];
    char output_folder[256];
    if ((argc != 3) && (argc !=2)) {
        printf ( "usage: make_demo2 input_folder output_folder.\
        \nNB input_folder must contain \"SpecificationFile.txt\", output will be wrtitten to \"output_folder/out_data_time/\".\
        \nIf output_folder is not given the value from SpecificationFile.txt will be used.\n" );
        return 0;
    } else {
        sprintf ( input_folder, "%s", argv[1] );
        sprintf ( output_folder, "%s", argv[2] );
        printf ( "input_folder = %s , output_folder = %s\n", input_folder, output_folder );
    }
    
    // Initialize
    //cuInit ( 0 );
    //int deviceCount = 0;
    size_t deviceCount;
    //CUDA: cuDeviceGetCount ( &deviceCount ); OpenCL:
    clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 0, NULL, &deviceCount);
    if ( deviceCount == 0 ) {
        printf ( "There is no device supporting CUDA.\n" );
        exit ( 0 );
    }
<<<<<<< HEAD:Morphogenesis/src/make_demo2.cpp
    cl_device_id devices;
    clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &devices, NULL);
    CUcontext clContext;
    clCreateContext ( &clContext, 0, devices );
=======
    cl_device_id clDevice;
    clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &clDevice, NULL);
    CUcontext cuContext;
    cuCtxCreate ( &cuContext, 0, clDevice );
>>>>>>> 75eada6585054e07bf9262a150f34af03aa68428:src/make_demo2.cpp
    
    FluidSystem fluid;
    fluid.InitializeOpenCL ();
    //std::cout<<"\n\nmake_demo2 chk0,"<<std::flush;
    
    fluid.ReadSpecificationFile ( input_folder );
    std::cout<<"\n\nmake_demo2 chk1, fluid.launchParams.debug="<<fluid.launchParams.debug<<", fluid.launchParams.genomePath=" <<fluid.launchParams.genomePath  << ",  fluid.launchParams.spacing="<<fluid.launchParams.spacing<<std::flush;
    
    for(int i=0; i<256; i++){fluid.launchParams.paramsPath[i] = input_folder[i];}
    for(int i=0; i<256; i++){fluid.launchParams.pointsPath[i] = input_folder[i];}
    //for(int i=0; i<256; i++){fluid.launchParams.genomePath[i] = input_folder[i];} // obtained from SpecificationFile.txt above.
    if(argc==3)for(int i=0; i<256; i++){fluid.launchParams.outPath[i] = output_folder[i];}
    
    if(mkdir(output_folder, 0755) == -1) cerr << "\nError :  failed to create output_folder.\n" << strerror(errno) << endl;
    else cout << "output_folder created\n"; // NB 0755 = rwx owner, rx for others.
    
    fluid.WriteDemoSimParams(           // Generates the simulation from data previously loaded from SpecificationFile.txt .
        fluid.launchParams.paramsPath, GPU_DUAL, CPU_YES, fluid.launchParams.num_particles, fluid.launchParams.spacing, fluid.launchParams.x_dim, fluid.launchParams.y_dim, fluid.launchParams.z_dim, fluid.launchParams.demoType, fluid.launchParams.simSpace, fluid.launchParams.debug
    ); /*const char * relativePath*/ 
    //std::cout<<"\n\nmake_demo2 chk2 "<<std::flush;
    uint num_particles_start=fluid.ActivePoints();
    
    fluid.TransferToCL (); 
    fluid.Run2Simulation ();
    
    //std::cout<<"\n\nmake_demo2 chk3 "<<std::flush;
    fluid.WriteResultsCSV(input_folder, output_folder, num_particles_start);// NB post-slurm script to (i) cat results.csv files, (ii)tar-gzip and ftp folders to recipient.
    
    size_t   free1, free2, total;
    cudaMemGetInfo(&free1, &total);
    printf("\n\nmake_demo2: Cuda Memory, before cuCtxDestroy(clContext): free=%lu, total=%lu.\t",free1,total);
    
<<<<<<< HEAD:Morphogenesis/src/make_demo2.cpp
    cl_int clResult = cuCtxDestroy ( clContext ) ;
=======
    cl_int clResult = cuCtxDestroy ( cuContext ) ;
>>>>>>> 75eada6585054e07bf9262a150f34af03aa68428:src/make_demo2.cpp
    if ( clResult!=0 ) {printf ( "error closing, clResult = %i \n",clResult );}
    
    cudaMemGetInfo(&free2, &total);
    printf("\nmake_demo2: After cuCtxDestroy(clContext): free=%lu, total=%lu, released=%lu.\n",free2,total,(free2-free1) );
    
    printf ( "\nClosed make_demo2.\n" );
    return 0;
}
