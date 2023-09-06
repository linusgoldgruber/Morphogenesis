// #ifndef RUNCL_H
// #define RUNCL_H
//
// #include <CL/cl.h>	//or <CL/opencl.hpp>
// #include <cstdio>
// #include <cstdlib>
// #include <cassert>
// #include <cstring>
// #include <iostream>
// #include <string>
// #include <fstream>
// #include <boost/filesystem.hpp>
// #include <boost/format.hpp>
// #include <jsoncpp/json/json.h>
// #include </home/goldi/opencv-4.x/modules/core/include/opencv2/core/hal/interface.h> //not the right way to include, otherwise I get an error though
// #include "fluid.h"
// 							// indices for float params passed to __const params_buf
// #define PIXELS			0	// TODO Can these be #included from a common header for both host and device code?
// #define ROWS			1
// #define COLS			2
// #define LAYERS			3
// #define MAX_INV_DEPTH	4
// #define MIN_INV_DEPTH	5
// #define INV_DEPTH_STEP	6
// #define ALPHA_G			7
// #define BETA_G			8	//  __kernel void CacheG4
// //#define EPSILON 		9	//  __kernel void UpdateQD		// epsilon = 0.1		//This is a conflicting definition with fluid_system.h:6
// #define SIGMA_Q 		10									// sigma_q = 0.0559017
// #define SIGMA_D 		11
// #define THETA			12
// #define LAMBDA			13	//  __kernel void UpdateA2
// #define SCALE_EAUX		14
//
// #define FUNC_INSERT			0
// 	#define	FUNC_COUNTING_SORT	1
// 	#define FUNC_QUERY			2
// 	#define FUNC_COMPUTE_PRESS	3
// 	#define FUNC_COMPUTE_FORCE	4
// 	#define FUNC_ADVANCE		5
// 	#define FUNC_EMIT			6
// 	#define FUNC_RANDOMIZE		7
// 	#define FUNC_SAMPLE			8
// 	#define FUNC_FPREFIXSUM		9
// 	#define FUNC_FPREFIXUP	10
// 	#define FUNC_TALLYLISTS     11
// 	#define FUNC_COMPUTE_DIFFUSION          12
// 	#define FUNC_COUNT_SORT_LISTS           13
// 	#define FUNC_COMPUTE_GENE_ACTION        14
// 	#define FUNC_TALLY_GENE_ACTION        35
// 	#define FUNC_COMPUTE_BOND_CHANGES       15
//
// 	#define FUNC_INSERT_CHANGES             16 //insertChanges
// 	#define FUNC_PREFIXUP_CHANGES           17 //prefixFixupChanges
// 	#define FUNC_PREFIXSUM_CHANGES          18 //prefixSumChanges
// 	#define FUNC_TALLYLISTS_CHANGES         19 //tally_changelist_lengths
// 	#define FUNC_COUNTING_SORT_CHANGES      20 //countingSortChanges
// 	#define FUNC_COMPUTE_NERVE_ACTION       21 //computeNerveActivation
//
// 	#define FUNC_COMPUTE_MUSCLE_CONTRACTION 22 //computeMuscleContraction
// 	#define FUNC_HEAL                       23 //heal
// 	#define FUNC_LENGTHEN_TISSUE            24 //lengthen_muscle
// 	#define FUNC_LENGTHEN_MUSCLE            25 //lengthen_tissue
// 	#define FUNC_SHORTEN_TISSUE             26 //shorten_muscle
// 	#define FUNC_SHORTEN_MUSCLE             27 //shorten_tissue
//
// 	#define FUNC_STRENGTHEN_TISSUE          28 //strengthen_muscle
// 	#define FUNC_STRENGTHEN_MUSCLE          29 //strengthen_tissue
// 	#define FUNC_WEAKEN_TISSUE              30 //weaken_muscle
// 	#define FUNC_WEAKEN_MUSCLE              31 //weaken_tissue
//
// 	#define FUNC_EXTERNAL_ACTUATION         32
// 	#define FUNC_FIXED                      33
// 	#define FUNC_CLEAN_BONDS                34
//
// 	#define FUNC_INIT_FCURAND_STATE         36
// 	#define FUNC_COUNTING_SORT_EPIGEN       37
//
// 	#define	FUNC_ASSEMBLE_MUSCLE_FIBRES_OUTGOING     38
// 	#define	FUNC_ASSEMBLE_MUSCLE_FIBRES_INCOMING     39
// 	#define	FUNC_INITIALIZE_BONDS                    40
// 	#define FUNC_MAX			            41
//
// using namespace std;
// class RunCL
// {
// public:
// 	Json::Value obj;
// 	int					verbosity;
// 	std::vector<cl_platform_id> m_platform_ids;
// 	cl_context			m_context;
// 	cl_device_id		m_device_id;
// 	cl_command_queue	m_queue, uload_queue, dload_queue, track_queue;
// 	cl_program			m_program;
// 	cl_kernel			m_Kern[FUNC_MAX];
// 	cl_mem				m_FParamDevice, m_FluidDevice, m_FluidTempDevice, m_FGenomeDevice;		//GPU pointer containers
//
// 	size_t  			global_work_size, local_work_size, image_size_bytes;
// 	bool 				gpu, amdPlatform;
// 	cl_device_id 		deviceId;
// 	float 				params[16] = {0};
// 	int 				width, height, costVolLayers, baseImage_type, count=0, keyFrameCount=0, costVolCount=0, QDcount=0, A_count=0;
// 	std::map< std::string, boost::filesystem::path > paths;
//
// 	RunCL(Json::Value obj_);
// 	void createFolders();	// Called by RunCL(..) constructor, above.
// 	void saveCostVols(float max_range);
// 	void updateQD(float epsilon, float theta, float sigma_q, float sigma_d);
// 	void updateA ( float lambda, float theta );
//
// 	void CleanUp();
// 	void allocatemem(FParams fparam, FBufs *fbuf, FBufs *ftemp, FGenome fgenome);
// 	void exit_(cl_int res);
// 	~RunCL();
//
// 	//void DownloadAndSaveVolume_3Channel(cl_mem buffer, std::string count, boost::filesystem::path folder, size_t image_size_bytes, cv::Size size_mat, int type_mat, bool show );
// 	//void initializeCostVol(float* k2k, cv::Mat &baseImage, cv::Mat &image, float *cdata, float *hdata, float thresh, int layers);
// 	//void initializeAD();
//
// 	int  convertToString(const char *filename, std::string& s){
// 		size_t size;
// 		char*  str;
// 		std::fstream f(filename, (std::fstream::in | std::fstream::binary));
// 		if (f.is_open() ) {
// 			size_t fileSize;
// 			f.seekg(0, std::fstream::end);
// 			size = fileSize = (size_t)f.tellg();
// 			f.seekg(0, std::fstream::beg);
// 			str = new char[size + 1];
// 			if (!str) {
// 				f.close();
// 				return 0;
// 			}
// 			f.read(str, fileSize);
// 			f.close();
// 			str[size] = '\0';
// 			s = str;
// 			delete[] str;
// 			return 0;
// 		}
// 											cout << "Error: failed to open file\n:" << filename << endl;
// 		return 1;
// 	}
//
// 	int waitForEventAndRelease(cl_event *event){
// 											if(verbosity>0) cout << "\nwaitForEventAndRelease_chk0, event="<<event<<" *event="<<*event << flush;
// 		cl_int status = CL_SUCCESS;
// 		status = clWaitForEvents(1, event); if (status != CL_SUCCESS) { cout << "\nclWaitForEvents status=" << status << ", " <<  checkerror(status) <<"\n" << flush; exit_(status); }
// 		status = clReleaseEvent(*event); 	if (status != CL_SUCCESS) { cout << "\nclReleaseEvent status="  << status << ", " <<  checkerror(status) <<"\n" << flush; exit_(status); }
// 		return status;
// 	}
//
// 	string checkerror(int input) {
// 		int errorCode = input;
// 		switch (errorCode) {
// 		case -9999:											return "Illegal read or write to a buffer";		// NVidia error code
// 		case CL_DEVICE_NOT_FOUND:							return "CL_DEVICE_NOT_FOUND";
// 		case CL_DEVICE_NOT_AVAILABLE:						return "CL_DEVICE_NOT_AVAILABLE";
// 		case CL_COMPILER_NOT_AVAILABLE:						return "CL_COMPILER_NOT_AVAILABLE";
// 		case CL_MEM_OBJECT_ALLOCATION_FAILURE:				return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
// 		case CL_OUT_OF_RESOURCES:							return "CL_OUT_OF_RESOURCES";
// 		case CL_OUT_OF_HOST_MEMORY:							return "CL_OUT_OF_HOST_MEMORY";
// 		case CL_PROFILING_INFO_NOT_AVAILABLE:				return "CL_PROFILING_INFO_NOT_AVAILABLE";
// 		case CL_MEM_COPY_OVERLAP:							return "CL_MEM_COPY_OVERLAP";
// 		case CL_IMAGE_FORMAT_MISMATCH:						return "CL_IMAGE_FORMAT_MISMATCH";
// 		case CL_IMAGE_FORMAT_NOT_SUPPORTED:					return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
// 		case CL_BUILD_PROGRAM_FAILURE:						return "CL_BUILD_PROGRAM_FAILURE";
// 		case CL_MAP_FAILURE:								return "CL_MAP_FAILURE";
// 		case CL_MISALIGNED_SUB_BUFFER_OFFSET:				return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
// 		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:	return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
// 		case CL_INVALID_VALUE:								return "CL_INVALID_VALUE";
// 		case CL_INVALID_DEVICE_TYPE:						return "CL_INVALID_DEVICE_TYPE";
// 		case CL_INVALID_PLATFORM:							return "CL_INVALID_PLATFORM";
// 		case CL_INVALID_DEVICE:								return "CL_INVALID_DEVICE";
// 		case CL_INVALID_CONTEXT:							return "CL_INVALID_CONTEXT";
// 		case CL_INVALID_QUEUE_PROPERTIES:					return "CL_INVALID_QUEUE_PROPERTIES";
// 		case CL_INVALID_COMMAND_QUEUE:						return "CL_INVALID_COMMAND_QUEUE";
// 		case CL_INVALID_HOST_PTR:							return "CL_INVALID_HOST_PTR";
// 		case CL_INVALID_MEM_OBJECT:							return "CL_INVALID_MEM_OBJECT";
// 		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:			return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
// 		case CL_INVALID_IMAGE_SIZE:							return "CL_INVALID_IMAGE_SIZE";
// 		case CL_INVALID_SAMPLER:							return "CL_INVALID_SAMPLER";
// 		case CL_INVALID_BINARY:								return "CL_INVALID_BINARY";
// 		case CL_INVALID_BUILD_OPTIONS:						return "CL_INVALID_BUILD_OPTIONS";
// 		case CL_INVALID_PROGRAM:							return "CL_INVALID_PROGRAM";
// 		case CL_INVALID_PROGRAM_EXECUTABLE:					return "CL_INVALID_PROGRAM_EXECUTABLE";
// 		case CL_INVALID_KERNEL_NAME:						return "CL_INVALID_KERNEL_NAME";
// 		case CL_INVALID_KERNEL_DEFINITION:					return "CL_INVALID_KERNEL_DEFINITION";
// 		case CL_INVALID_KERNEL:								return "CL_INVALID_KERNEL";
// 		case CL_INVALID_ARG_INDEX:							return "CL_INVALID_ARG_INDEX";
// 		case CL_INVALID_ARG_VALUE:							return "CL_INVALID_ARG_VALUE";
// 		case CL_INVALID_ARG_SIZE:							return "CL_INVALID_ARG_SIZE";
// 		case CL_INVALID_KERNEL_ARGS:						return "CL_INVALID_KERNEL_ARGS";
// 		case CL_INVALID_WORK_DIMENSION:						return "CL_INVALID_WORK_DIMENSION";
// 		case CL_INVALID_WORK_GROUP_SIZE:					return "CL_INVALID_WORK_GROUP_SIZE";
// 		case CL_INVALID_WORK_ITEM_SIZE:						return "CL_INVALID_WORK_ITEM_SIZE";
// 		case CL_INVALID_GLOBAL_OFFSET:						return "CL_INVALID_GLOBAL_OFFSET";
// 		case CL_INVALID_EVENT_WAIT_LIST:					return "CL_INVALID_EVENT_WAIT_LIST";
// 		case CL_INVALID_EVENT:								return "CL_INVALID_EVENT";
// 		case CL_INVALID_OPERATION:							return "CL_INVALID_OPERATION";
// 		case CL_INVALID_GL_OBJECT:							return "CL_INVALID_GL_OBJECT";
// 		case CL_INVALID_BUFFER_SIZE:						return "CL_INVALID_BUFFER_SIZE";
// 		case CL_INVALID_MIP_LEVEL:							return "CL_INVALID_MIP_LEVEL";
// 		case CL_INVALID_GLOBAL_WORK_SIZE:					return "CL_INVALID_GLOBAL_WORK_SIZE";
// #if CL_HPP_MINIMUM_OPENCL_VERSION >= 200
// 		case CL_INVALID_DEVICE_QUEUE:						return "CL_INVALID_DEVICE_QUEUE";
// 		case CL_INVALID_PIPE_SIZE:							return "CL_INVALID_PIPE_SIZE";
// #endif
// 		default:											return "unknown error code";
// 		}
// 	}
//
// 	void InitializeKernels(cl_program m_program) {
// 	// Setup the array of kernels
// 		m_Kern[FUNC_INSERT] = clCreateKernel(m_program, "insertParticles", NULL);
// 		m_Kern[FUNC_COUNTING_SORT] = clCreateKernel(m_program, "countingSortFull", NULL);
// 		m_Kern[FUNC_QUERY] = clCreateKernel(m_program, "computeQuery", NULL);
// 		m_Kern[FUNC_COMPUTE_PRESS] = clCreateKernel(m_program, "computePressure", NULL);
// 		m_Kern[FUNC_COMPUTE_FORCE] = clCreateKernel(m_program, "computeForce", NULL);
// 		m_Kern[FUNC_ADVANCE] = clCreateKernel(m_program, "advanceParticles", NULL);
// 		m_Kern[FUNC_EMIT] = clCreateKernel(m_program, "emitParticles", NULL);
// 		m_Kern[FUNC_RANDOMIZE] = clCreateKernel(m_program, "randomInit", NULL);
// 		m_Kern[FUNC_SAMPLE] = clCreateKernel(m_program, "sampleParticles", NULL);
// 		m_Kern[FUNC_FPREFIXSUM] = clCreateKernel(m_program, "prefixSum", NULL);
// 		m_Kern[FUNC_FPREFIXUP] = clCreateKernel(m_program, "prefixFixup", NULL);
// 		m_Kern[FUNC_TALLYLISTS] = clCreateKernel(m_program, "tally_denselist_lengths", NULL);
// 		m_Kern[FUNC_COMPUTE_DIFFUSION] = clCreateKernel(m_program, "computeDiffusion", NULL);
// 		m_Kern[FUNC_COUNT_SORT_LISTS] = clCreateKernel(m_program, "countingSortDenseLists", NULL);
// 		m_Kern[FUNC_COMPUTE_GENE_ACTION] = clCreateKernel(m_program, "computeGeneAction", NULL);
// 		m_Kern[FUNC_TALLY_GENE_ACTION] = clCreateKernel(m_program, "tallyGeneAction", NULL);
// 		m_Kern[FUNC_COMPUTE_BOND_CHANGES] = clCreateKernel(m_program, "computeBondChanges", NULL);
// 		m_Kern[FUNC_COUNTING_SORT_CHANGES] = clCreateKernel(m_program, "countingSortChanges", NULL);
// 		m_Kern[FUNC_COMPUTE_NERVE_ACTION] = clCreateKernel(m_program, "computeNerveActivation", NULL);
// 		m_Kern[FUNC_COMPUTE_MUSCLE_CONTRACTION] = clCreateKernel(m_program, "computeMuscleContraction", NULL);
// 		m_Kern[FUNC_CLEAN_BONDS] = clCreateKernel(m_program, "cleanBonds", NULL);
// 		m_Kern[FUNC_HEAL] = clCreateKernel(m_program, "heal", NULL);
// 		m_Kern[FUNC_LENGTHEN_MUSCLE] = clCreateKernel(m_program, "lengthen_muscle", NULL);
// 		m_Kern[FUNC_LENGTHEN_TISSUE] = clCreateKernel(m_program, "lengthen_tissue", NULL);
// 		m_Kern[FUNC_SHORTEN_MUSCLE] = clCreateKernel(m_program, "shorten_muscle", NULL);
// 		m_Kern[FUNC_SHORTEN_TISSUE] = clCreateKernel(m_program, "shorten_tissue", NULL);
// 		m_Kern[FUNC_STRENGTHEN_MUSCLE] = clCreateKernel(m_program, "strengthen_muscle", NULL);
// 		m_Kern[FUNC_STRENGTHEN_TISSUE] = clCreateKernel(m_program, "strengthen_tissue", NULL);
// 		m_Kern[FUNC_WEAKEN_MUSCLE] = clCreateKernel(m_program, "weaken_muscle", NULL);
// 		m_Kern[FUNC_WEAKEN_TISSUE] = clCreateKernel(m_program, "weaken_tissue", NULL);
// 		m_Kern[FUNC_EXTERNAL_ACTUATION] = clCreateKernel(m_program, "externalActuation", NULL);
// 		m_Kern[FUNC_FIXED] = clCreateKernel(m_program, "fixedParticles", NULL);
// 		m_Kern[FUNC_INIT_FCURAND_STATE] = clCreateKernel(m_program, "initialize_FCURAND_STATE", NULL);
// 		m_Kern[FUNC_ASSEMBLE_MUSCLE_FIBRES_OUTGOING] = clCreateKernel(m_program, "assembleMuscleFibresOutGoing", NULL);
// 		m_Kern[FUNC_ASSEMBLE_MUSCLE_FIBRES_INCOMING] = clCreateKernel(m_program, "assembleMuscleFibresInComing", NULL);
// 		m_Kern[FUNC_INITIALIZE_BONDS] = clCreateKernel(m_program, "initialize_bonds", NULL);
// 	}
// 	/*
// 	void ReadOutput(uchar* outmat) {
// 		ReadOutput(outmat, amem,  (width * height * sizeof(float)) );
// 	}
//
// 	void ReadOutput(uchar* outmat, cl_mem buf_mem, size_t data_size, size_t offset=0) {
// 		cl_event readEvt;
// 		cl_int status;
// 														cout<<"\nReadOutput: &outmat="<<&outmat<<", buf_mem="<<buf_mem<<", data_size="<<data_size<<", offset="<<offset<<"\t"<<flush;
// 		status = clEnqueueReadBuffer(dload_queue,			// command_queue
// 											buf_mem,		// buffer
// 											CL_FALSE,		// blocking_read
// 											offset,			// offset
// 											data_size,		// size
// 											outmat,			// pointer
// 											0,				// num_events_in_wait_list
// 											NULL,			// event_waitlist
// 											&readEvt);		// event
// 														if (status != CL_SUCCESS) { cout << "\nclEnqueueReadBuffer(..) status=" << checkerror(status) <<"\n"<<flush; exit_(status);}
// 															else if(verbosity>0) cout <<"\nclEnqueueReadBuffer(..)"<<flush;
// 		status = clFlush(dload_queue);					if (status != CL_SUCCESS) { cout << "\nclFlush(m_queue) status = " 		<< checkerror(status) <<"\n"<<flush; exit_(status);}
// 															else if(verbosity>0) cout <<"\nclFlush(..)"<<flush;
// 		status = clWaitForEvents(1, &readEvt); 			if (status != CL_SUCCESS) { cout << "\nclWaitForEvents status="			<< checkerror(status) <<"\n"<<flush; exit_(status);}
// 															else if(verbosity>0) cout <<"\nclWaitForEvents(..)"<<flush;
// 	}*/
// };
//
// #endif /*RUNCL_H*/
