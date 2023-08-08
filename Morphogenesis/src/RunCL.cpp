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
	createFolders( );																	/*Step1: Getting platforms and choose an available one.*/////////
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
																						// Multiple queues for latency hiding: Upload, Download, Mapping, Tracking,... autocalibration, SIRFS, SPMP
																						// NB Might want to create command queues on multiple platforms & devices.
																						// NB might want to divde a task across multiple MPI Ranks on a multi-GPU WS or cluster.
	const char *filename = obj["kernel_filepath"].asCString();							/*Step 5: Create program object*///////////////////////////////
	string sourceStr;
	status 						= convertToString(filename, sourceStr);					if(status!=CL_SUCCESS)	{cout<<"\nstatus="<<checkerror(status)<<"\n"<<flush;exit_(status);}
	const char 	*source 		= sourceStr.c_str();
	size_t 		sourceSize[] 	= { strlen(source) };
	m_program 	= clCreateProgramWithSource(m_context, 1, &source, sourceSize, NULL);
	
	status = clBuildProgram(m_program, 1, devices, NULL, NULL, NULL);					/*Step 6: Build program.*/////////////////////
	if (status != CL_SUCCESS){
		printf("\nclBuildProgram failed: %d\n", status);
		char buf[0x10000];
		clGetProgramBuildInfo(m_program, deviceId, CL_PROGRAM_BUILD_LOG, 0x10000, buf, NULL);
		printf("\n%s\n", buf);
		exit_(status);
	}
	cost_kernel     = clCreateKernel(m_program, "BuildCostVolume2", NULL);				/*Step 7: Create kernel objects.*////////////
	cache3_kernel   = clCreateKernel(m_program, "CacheG3", 			NULL);
	updateQD_kernel = clCreateKernel(m_program, "UpdateQD", 		NULL);
	updateA_kernel  = clCreateKernel(m_program, "UpdateA2", 		NULL);
	basemem=imgmem=cdatabuf=hdatabuf=k2kbuf=dmem=amem=basegraymem=gxmem=gymem=g1mem=lomem=himem=0;		// set device pointers to zero
																						if(verbosity>0) cout << "RunCL_constructor finished\n" << flush;
}

void RunCL::createFolders(){
																						if(verbosity>0) cout << "\n createFolders_chk 0\n" << flush;
	std::time_t   result  = std::time(nullptr);
	std::string   out_dir = std::asctime(std::localtime(&result));
	out_dir.pop_back(); 																// req to remove new_line from end of string.

	boost::filesystem::path 	out_path(boost::filesystem::current_path());
	boost::filesystem::path 	conf_outpath( obj["outpath"].asString() );
	if (conf_outpath.empty() ||  conf_outpath.is_absolute() ) {
		out_path = out_path.parent_path().parent_path();								// move "out_path" up two levels in the directory tree.
		out_path += conf_outpath;
	}else {out_path = conf_outpath;}
	out_path += "/output/";
	if(boost::filesystem::create_directory(out_path)) { std::cerr<< "Directory Created: "<<out_path<<std::endl;}else{ std::cerr<< "Output directory previously created: "<<out_path<<std::endl;}
	out_path +=  out_dir;																if(verbosity>0) cout <<"Creating output sub-directories: "<< out_path <<std::endl;
	boost::filesystem::create_directory(out_path);
	out_path += "/";																	if(verbosity>0) cout << "\n createFolders_chk 1\n" << flush;

	boost::filesystem::path temp_path = out_path;								// Vector of device buffer names
	std::vector<std::string> names = {"basemem","imgmem","cdatabuf","hdatabuf","pbuf","dmem", "amem","basegraymem","gxmem","gymem","g1mem","qmem","lomem","himem", "img_sum_buf"};
	std::pair<std::string, boost::filesystem::path> tempPair;

	for (std::string key : names){
		temp_path = out_path;
		temp_path += key;
		tempPair = {key, temp_path};
		paths.insert(tempPair);
		boost::filesystem::create_directory(temp_path);
		temp_path += "/png/";
		boost::filesystem::create_directory(temp_path);
	}
	if(verbosity>0) {
		cout << "\nRunCL::createFolders() chk1\n";
		cout << "KEY\tPATH\n";														// print the folder paths
		for (auto itr=paths.begin(); itr!=paths.end(); ++itr) { cout<<"First:["<< itr->first << "]\t:\t Second:"<<itr->second<<"\n"; }
		cout<<"\npaths.at(\"basemem\")="<<paths.at("basemem")<<"\n"<<flush;
	}
}

void RunCL::saveCostVols(float max_range)
{																				if(verbosity>0) cout<<"\nsaveCostVols: Calling DownloadAndSaveVolume";
	stringstream ss;
	ss << "saveCostVols";
	ss << (keyFrameCount*1000 + costVolCount);
	DownloadAndSaveVolume(cdatabuf, 	ss.str(), paths.at("cdatabuf"), 	width * height * sizeof(float), baseImage_size, CV_32FC1,  false  , max_range);
	DownloadAndSaveVolume(hdatabuf, 	ss.str(), paths.at("hdatabuf"), 	width * height * sizeof(float), baseImage_size, CV_32FC1,  false  , max_range);
	DownloadAndSaveVolume(img_sum_buf, 	ss.str(), paths.at("img_sum_buf"), 	width * height * sizeof(float), baseImage_size, CV_32FC1,  false  , max_range);
																				if(verbosity>0) cout <<"\ncostVolCount="<<costVolCount << "\ncalcCostVol chk13_finished\n" << flush;
}

/*
buffers being saved:

DownloadAndSave:
amem CV_32FC1	max 1
qmem			max 1
dmem			max 1
lomem			max 1
himem			max 1
basegraymem     max 1
gxmem			max 1
gymem			max 1
g1mem			max 1

imgmem CV_32FC3 max(1,1,1)
basemem			max(1,1,1)

DownloadAndSaveVolume:
cdatabuf CV_32FC1		max 1
hdatabuf				max 30.001
 */
void RunCL::DownloadAndSave(cl_mem buffer, std::string count, boost::filesystem::path folder_tiff, size_t image_size_bytes, cv::Size size_mat, int type_mat, bool show, float max_range ){
																				if(verbosity>0) cout<<"\n\nDownloadAndSave filename = ["<<folder_tiff.filename().string()<<"] folder="<<folder_tiff<<", image_size_bytes="<<image_size_bytes<<", size_mat="<<size_mat<<", type_mat="<<size_mat<<"\t"<<flush;
		cv::Mat temp_mat = cv::Mat::zeros (size_mat, type_mat);					// (int rows, int cols, int type)
		ReadOutput(temp_mat.data, buffer,  image_size_bytes); 					// NB contains elements of type_mat, (CV_32FC1 for most buffers)
		cv::Scalar sum = cv::sum(temp_mat);										// NB always returns a 4 element vector.

		double minVal=1, maxVal=1;
		cv::Point minLoc={0,0}, maxLoc{0,0};
		if (temp_mat.channels()==1) { cv::minMaxLoc(temp_mat, &minVal, &maxVal, &minLoc, &maxLoc); }
		string type_string = checkCVtype(type_mat);
		stringstream ss;
		stringstream png_ss;
		ss << "/" << folder_tiff.filename().string() << "_" << count <<"_sum"<<sum<<"type_"<<type_string<<"min"<<minVal<<"max"<<maxVal<<"maxRange"<<max_range;
		png_ss << "/" << folder_tiff.filename().string() << "_" << count;
		boost::filesystem::path folder_png = folder_tiff;
		folder_tiff += ss.str();
		folder_tiff += ".tiff";
		folder_png  += "/png/";

		folder_png  += png_ss.str();
		folder_png  += ".png";
																				if(verbosity>0) cout<<"\n\nDownloadAndSave filename = ["<<ss.str()<<"]";
		cv::Mat outMat;
		if (type_mat != CV_32FC1) {
			cout << "Error  (type_mat != CV_32FC1)" << flush;
			return;
		}
		if (max_range == 0){ temp_mat /= maxVal;}								// Squash/stretch & shift to 0.0-1.0 range
		else if (max_range <0.0){
			temp_mat /=(-2*max_range);
			temp_mat +=0.5;
		}else{ temp_mat /=max_range;}

		cv::imwrite(folder_tiff.string(), temp_mat );
		temp_mat *= 256*256;
		temp_mat.convertTo(outMat, CV_16UC1);
		cv::imwrite(folder_png.string(), outMat );
		if(show) cv::imshow( ss.str(), outMat );
}

void RunCL::DownloadAndSave_3Channel(cl_mem buffer, std::string count, boost::filesystem::path folder_tiff, size_t image_size_bytes, cv::Size size_mat, int type_mat, bool show ){
																				if(verbosity>0) cout<<"\n\nDownloadAndSave_3Channel filename = ["<<folder_tiff.filename()<<"] folder="<<folder_tiff<<", image_size_bytes="<<image_size_bytes<<", size_mat="<<size_mat<<", type_mat="<<type_mat<<"\t"<<flush;
		cv::Mat temp_mat = cv::Mat::zeros (size_mat, type_mat);
		ReadOutput(temp_mat.data, buffer,  image_size_bytes);
		cv::Scalar 	sum = cv::sum(temp_mat);									// NB always returns a 4 element vector.
		string 		type_string=checkCVtype(type_mat);
		double 		minVal[3]={1,1,1}, 					maxVal[3]={0,0,0};
		cv::Point 	minLoc[3]={{0,0},{0,0},{0,0}}, 		maxLoc[3]={{0,0},{0,0},{0,0}};
		vector<cv::Mat> spl;
		split(temp_mat, spl);													// process - extract only the correct channel
		double max = 0;
		for (int i =0; i < 3; ++i){
			cv::minMaxLoc(spl[i], &minVal[i], &maxVal[i], &minLoc[i], &maxLoc[i]);
			if (maxVal[i] > max) max = maxVal[i];
		}
		stringstream ss;
		stringstream png_ss;
		ss<<"/"<<folder_tiff.filename().string()<<"_"<<count<<"_sum"<<sum<<"type_"<<type_string<<"min("<<minVal[0]<<","<<minVal[1]<<","<<minVal[2]<<")_max("<<maxVal[0]<<","<<maxVal[1]<<","<<maxVal[2]<<")";
		png_ss<< "/" << folder_tiff.filename().string() << "_" << count;
		if(show){
			cv::Mat temp;
			temp_mat.convertTo(temp, CV_8U);									// NB need CV_U8 for imshow(..)
			cv::imshow( ss.str(), temp);
		}
		boost::filesystem::path folder_png = folder_tiff;
		folder_png  += "/png/";
		folder_png  += png_ss.str();
		folder_png  += ".png";

		folder_tiff += ss.str();
		folder_tiff += ".tiff";

		cv::Mat outMat;
		if (type_mat == CV_32FC3){
				cv::imwrite(folder_tiff.string(), temp_mat );
				temp_mat *=256;
				temp_mat.convertTo(outMat, CV_8U);
				cv::imwrite(folder_png.string(), (outMat) );					// Has "Grayscale 16-bit gamma integer"
		}
}

void RunCL::DownloadAndSaveVolume(cl_mem buffer, std::string count, boost::filesystem::path folder, size_t image_size_bytes, cv::Size size_mat, int type_mat, bool show, float max_range ){
																				if(verbosity>0) {
																					cout<<"\n\nDownloadAndSaveVolume, costVolLayers="<<costVolLayers<<", filename = ["<<folder.filename().string()<<"]";
																					cout<<"\n folder="<<folder.string()<<",\t image_size_bytes="<<image_size_bytes<<",\t size_mat="<<size_mat<<",\t type_mat="<<size_mat<<"\t"<<flush;
																				}
	cv::Mat temp_mat = cv::Mat::zeros (size_mat, type_mat);						//(int rows, int cols, int type)

	for(int i=0; i<costVolLayers; i++){
																				if(verbosity>0) cout << "\ncostVolLayers="<<costVolLayers<<", i="<<i<<"\t";
		size_t offset = i * image_size_bytes;
		ReadOutput(temp_mat.data, buffer,  image_size_bytes, offset);

		cv::Scalar sum = cv::sum(temp_mat);
		double minVal=1, maxVal=1;
		cv::Point minLoc={0,0}, maxLoc{0,0};
		if (type_mat == CV_32FC1) cv::minMaxLoc(temp_mat, &minVal, &maxVal, &minLoc, &maxLoc);

		boost::filesystem::path new_filepath = folder;
		boost::filesystem::path folder_png   = folder;

		string type_string = checkCVtype(type_mat);
		stringstream ss;
		stringstream png_ss;
		ss << "/"<< folder.filename().string() << "_" << count << "_layer"<< i <<"_sum"<<sum<<"type_"<<type_string<< "min"<<minVal<<"max"<<maxVal;
		png_ss << "/"<< folder.filename().string() << "_" << count << "_layer"<< i;
		if(show){
			cv::Mat temp;
			temp_mat.convertTo(temp, CV_8U);									// NB need CV_U8 for imshow(..)
			cv::imshow(ss.str(), temp);
		}
		new_filepath += ss.str();
		new_filepath += ".tiff";
		folder_png += "/png/";
		folder_png += png_ss.str();
		folder_png += ".png";
																				if(verbosity>0) cout << "\nnew_filepath.string() = "<<new_filepath.string() <<"\n";
		cv::Mat outMat;

		if (type_mat != CV_32FC1) {
			cout << "Error  (type_mat != CV_32FC1)" << flush;
			return;
		}
		if (max_range == 0){ temp_mat /= maxVal;}								// Squash/stretch & shift to 0.0-1.0 range
		else if (max_range <0.0){
			temp_mat /=(-2*max_range);
			temp_mat +=0.5;
		}else{ temp_mat /=max_range;}

		cv::imwrite(new_filepath.string(), temp_mat );
		temp_mat *= 256*256;
		temp_mat.convertTo(outMat, CV_16UC1);
		cv::imwrite(folder_png.string(), outMat );
		if(show) cv::imshow( ss.str(), outMat );
	}
}

void RunCL::allocatemem(float* gx, float* gy, float* params, int layers, cv::Mat &baseImage, float *cdata, float *hdata, float *img_sum_data)
{
																				if(verbosity>0) cout << "RunCL::allocatemem_chk0" << flush;
	stringstream 		ss;
	ss << "allocatemem";
	cl_int status;
	cl_event writeEvt;

	image_size_bytes	= baseImage.total() * baseImage.elemSize() ;
	costVolLayers 		= layers;
	baseImage_size 		= baseImage.size();
	baseImage_type 		= baseImage.type();
	int layerstep 		= width * height;
	// Get the maximum work group size for executing the kernel on the device ///////// From https://github.com/rsnemmen/OpenCL-examples/blob/e2c34f1dfefbd265cfb607c2dd6c82c799eb322a/square_array/square.c
	status = clGetKernelWorkGroupInfo(cost_kernel, deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local_work_size), &local_work_size, NULL); 	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}

	// Number of total work items, calculated here after 1st image is loaded &=> know the size.
	// NB localSize must be devisor
	// NB global_work_size must be a whole number of "Preferred work group size multiple" for Nvidia.
	global_work_size = ceil((float)layerstep/(float)local_work_size) * local_work_size;
	local_work_size=32; // trial for nvidia
																				if(verbosity>1){
																					cout<<"\nglobal_work_size="<<global_work_size<<", local_work_size="<<local_work_size<<", deviceId="<<deviceId<<"\n"<<flush;
																					cout<<"\nlayerstep=width*height="<<width<<"*"<<height<<"="<<layerstep<<",\tsizeof(layerstep)="<< sizeof(layerstep) <<",\tsizeof(int)="<< sizeof(int) <<flush;
																					cout<<"\n";
																					cout<<"\nallocatemem chk1, baseImage.total()=" << baseImage.total() << ", sizeof(float)="<< sizeof(float)<<flush;
																					cout<<"\nbaseImage.elemSize()="<< baseImage.elemSize()<<", baseImage.elemSize1()="<<baseImage.elemSize1()<<flush;
																					cout<<"\nbaseImage.type()="<< baseImage.type() <<", sizeof(baseImage.type())="<< sizeof(baseImage.type())<<flush;
																					cout<<"\n";
																					cout<<"\nallocatemem chk2, image_size_bytes="<< image_size_bytes <<  ", sizeof(float)="<< sizeof(float)<<flush;
																				}
	cl_int res;
	dmem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	amem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	gxmem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	gymem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	qmem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , 2*width*height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	g1mem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	lomem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	himem		= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	param_buf	= clCreateBuffer(m_context, CL_MEM_READ_ONLY  , 		    16 * sizeof(float), 0, &res);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}

	cdatabuf	= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * layers * sizeof(float), 0, &res);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	hdatabuf 	= clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * layers * sizeof(float), 0, &res);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	imgmem		= clCreateBuffer(m_context, CL_MEM_READ_ONLY  , 					   image_size_bytes, 0, &res);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	basemem		= clCreateBuffer(m_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,  image_size_bytes, 0, &res);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	k2kbuf		= clCreateBuffer(m_context, CL_MEM_READ_ONLY  ,  					   16*sizeof(float), 0, &res); 		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	img_sum_buf = clCreateBuffer(m_context, CL_MEM_READ_WRITE , width * height * layers * sizeof(float), 0, &res);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}

																				if(verbosity>1) {
																					cout << ",dmem = " 		<< dmem << endl;
																					cout << ",amem = " 		<< amem << endl;
																					cout << ",gxmem = " 	<< gxmem << endl;
																					cout << ",gymem = " 	<< gymem << endl;
																					cout << ",qmem = " 		<< qmem << endl;
																					cout << ",g1mem = " 	<< g1mem << endl;
																					cout << ",lomem = " 	<< lomem << endl;
																					cout << ",himem = " 	<< himem << endl;
																					cout << ",param_buf = " << param_buf << endl;
																					cout << ",cdatabuf = " 	<< cdatabuf << endl;
																					cout << ",hdatabuf = " 	<< hdatabuf << endl;
																					cout << ",imgmem = " 	<< imgmem << endl;
																					cout << ",basemem = " 	<< basemem << endl;
																					cout << ",k2kbuf = " 	<< k2kbuf << endl;
																					cout << ",imgmem = "	<< imgmem <<flush;
																					cout << ",imgmem image_size_bytes= "<< image_size_bytes <<flush;
																				}
	cv::Mat img_sum;
	baseImage.convertTo(img_sum, CV_32FC1);

	status = clEnqueueWriteBuffer(uload_queue, gxmem, 		CL_FALSE, 0, width*height*sizeof(float), 		gx, 			0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.3\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, gymem, 		CL_FALSE, 0, width*height*sizeof(float), 		gy, 			0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.4\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, param_buf, 	CL_FALSE, 0, 16 * sizeof(float), 				params, 		0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.5\n" << endl; exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, basemem, 	CL_FALSE, 0, width*height*sizeof(float)*3, 		baseImage.data, 0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.6\n" << endl;exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, cdatabuf, 	CL_FALSE, 0, width*height*layers*sizeof(float), cdata, 			0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.8\n" << endl;exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, hdatabuf, 	CL_FALSE, 0, width*height*layers*sizeof(float), hdata, 			0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.9\n" << endl;exit_(status);}
	status = clEnqueueWriteBuffer(uload_queue, img_sum_buf, CL_FALSE, 0, width*height*layers*sizeof(float), img_sum_data, 	0, NULL, &writeEvt);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: allocatemem_chk1.10\n" << endl;exit_(status);}
	clFlush(uload_queue); status = clFinish(uload_queue); 					if (status != CL_SUCCESS)	{ cout << "\nclFinish(uload_queue)="				<< status << checkerror(status)<<"\n"<<flush; exit_(status);}

	DownloadAndSave_3Channel(basemem, ss.str(), paths.at("basemem"), image_size_bytes, baseImage_size, baseImage_type, false );				// DownloadAndSave_3Channel(basemem,..) verify uploads.
																																			// set kernelArg. NB "0 &k2kbuf" & "2 &imgmem" set in calcCostVol(..)
	res = clSetKernelArg(cost_kernel, 1, sizeof(cl_mem),  &basemem);		if(res!=CL_SUCCESS){cout<<"\nbasemem res= "   		<<checkerror(res)<<"\n"<<flush;exit_(res);} // base
	res = clSetKernelArg(cost_kernel, 3, sizeof(cl_mem),  &cdatabuf);		if(res!=CL_SUCCESS){cout<<"\ncdatabuf res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} // cdata
	res = clSetKernelArg(cost_kernel, 4, sizeof(cl_mem),  &hdatabuf);		if(res!=CL_SUCCESS){cout<<"\nhdatabuf res = " 		<<checkerror(res)<<"\n"<<flush;exit_(res);} // hdata
	res = clSetKernelArg(cost_kernel, 5, sizeof(cl_mem),  &lomem);			if(res!=CL_SUCCESS){cout<<"\nlomem res = "    		<<checkerror(res)<<"\n"<<flush;exit_(res);} // lo
	res = clSetKernelArg(cost_kernel, 6, sizeof(cl_mem),  &himem);			if(res!=CL_SUCCESS){cout<<"\nhimem res = "    		<<checkerror(res)<<"\n"<<flush;exit_(res);} // hi
	res = clSetKernelArg(cost_kernel, 7, sizeof(cl_mem),  &amem);			if(res!=CL_SUCCESS){cout<<"\namem res = "     		<<checkerror(res)<<"\n"<<flush;exit_(res);} // a
	res = clSetKernelArg(cost_kernel, 8, sizeof(cl_mem),  &dmem);			if(res!=CL_SUCCESS){cout<<"\ndmem res = "     		<<checkerror(res)<<"\n"<<flush;exit_(res);} // d
	res = clSetKernelArg(cost_kernel, 9, sizeof(cl_mem),  &param_buf);		if(res!=CL_SUCCESS){cout<<"\nparam_buf res = "		<<checkerror(res)<<"\n"<<flush;exit_(res);} // param_buf
	res = clSetKernelArg(cost_kernel,10, sizeof(cl_mem),  &img_sum_buf);	if(res!=CL_SUCCESS){cout<<"\nimg_sum_buf res = " 	<<checkerror(res)<<"\n"<<flush;exit_(res);} // cdata
																				if(verbosity>0) cout << "RunCL::allocatemem_finished\n\n" << flush;
}

void RunCL::calcCostVol(float* k2k,  cv::Mat &image)
{
																				if(verbosity>0) cout << "\ncalcCostVol chk0," << flush;
	if (basemem == 0) { cout<<"\n\nERROR: calcCostVol basemem _not_ yet allocated\n\n"<<flush; exit(0); }
	stringstream ss;
	ss << "calcCostVol";
	size_t old_image_size_bytes = image_size_bytes;
	image_size_bytes 			= image.total() * image.elemSize() ;																		// num array elements * bytes per element all channels
	costVolCount++;
	keyFrameCount++;
	if (costVolCount > 1 && old_image_size_bytes!=image_size_bytes){cout << "\n\nMismatch 'image_size_bytes' !!!\n"<<flush; exit_(0); }
	cl_int status, res;
	cl_event writeEvt, ev;
																				if(verbosity>1) {
																					cout <<"\ncalcCostVol chk1, sizeof(float)="<< sizeof(float)<<"\n"<<flush;
																					cout <<"\nimage_size_bytes="  << image_size_bytes <<", image.total()="       <<image.total()        <<", sizeof(float)="<< sizeof(float)<<flush;
																					cout <<"\nimage.elemSize()="  << image.elemSize() <<", image.elemSize1()="   <<image.elemSize1()    <<flush;
																					cout <<"\nImage.type()="      << image.type()     <<", sizeof(Image.type())="<< sizeof(image.type())<<flush;
																					cout <<"\nm_queue=" 		  << m_queue 		  << "\nuload_queue=" 		 <<uload_queue;
																					cout <<"\nCL_FALSE="		  <<CL_FALSE          << "\n&writeEvt=" 		 <<&writeEvt			<<flush;
																					cout <<"\nimgmem = "		  <<imgmem            <<flush;
																					cout <<"\nimgmem image_size_bytes= "<< image_size_bytes <<flush;
																					cout <<"\nimage.size=="<<image.size <<", width="<<width <<", height="<<height <<flush;
																				}
	status = clEnqueueWriteBuffer(uload_queue, imgmem, CL_FALSE, 0, image_size_bytes, image.data, 0, NULL, &writeEvt);							// WriteBuffer imgmem #########
	if (status != CL_SUCCESS)	{ cout << "\nclEnqueueWriteBuffer imgmem status = " << checkerror(status) <<"\n"<<flush; exit_(status);}

	status = clEnqueueWriteBuffer(uload_queue, k2kbuf, CL_FALSE, 0, 16*sizeof(float), k2k,  0, NULL, &writeEvt);								// WriteBuffer k2k #### cam2cam pixel transform
	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}

	status = clFlush(uload_queue); 				if (status != CL_SUCCESS)	{ cout << "\nclFlush(m_queue) status = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clFinish(uload_queue); 			if (status != CL_SUCCESS)	{ cout << "\nclFinish(m_queue)="<<status<<" "<<checkerror(status)<<"\n"<<flush; exit_(status);}

																				if(verbosity>1) cout << "\ncalcCostVol chk10," << flush;		// set kernelArg
	res = clSetKernelArg(cost_kernel, 0, sizeof(cl_mem), &k2kbuf);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(cost_kernel, 2, sizeof(cl_mem), &imgmem);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
																				if(verbosity>1) {
																					cout<<"\n\npose transformation matrix, as floats:\n";
																					for(int p_iter=0;p_iter<16;p_iter++)  { cout<<"rt["<<p_iter<<"]="<<k2k[p_iter]<<"\t\t";   if((p_iter+1)%4==0){cout<<"\n"; };  }
																					cout<<"\n";
																					cout<<"\nk2kbuf="<<k2kbuf;
																					cout<<"\n\nclEnqueueNDRang*eKernel(";
																					cout<<"\ncommand_queue="<<m_queue;
																					cout<<"\nkernel=cost_kernel="<<cost_kernel;
																					cout<<"\nwork_dim="<<1;
																					cout<<"\nglobal_work_offset="<<"NULL";
																					cout<<"\nglobal_work_size="<<global_work_size;
																					cout<<"\nlocal_work_size="<<local_work_size;
																					cout<<"\nnum_events_in_wait_list"<<0;
																					cout<<"\nevent_wait_list=NULL";
																					cout<<"\nevent="<<ev;
																					cout<<"\n)"<<flush;
																				}
	status = clFlush(m_queue); 				if (status != CL_SUCCESS)	{ cout << "\nclFlush(m_queue) status = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clFinish(m_queue); 			if (status != CL_SUCCESS)	{ cout << "\nclFinish(m_queue)="		 <<checkerror(status)  <<"\n"<<flush; exit_(status);}

	res    = clEnqueueNDRangeKernel(m_queue, cost_kernel, 1, 0, &global_work_size, 0, 0, NULL, &ev); 										// ####### cost_kernel ############
	if (res != CL_SUCCESS)	{ cout << "\nclEnqueueNDRangeKernel res = " << checkerror(res) <<"\n"<<flush; exit_(res);}
	status = clFlush(m_queue); 				if (status != CL_SUCCESS)	{ cout << "\nclFlush(m_queue) status = " <<checkerror(status)  <<"\n"<<flush; exit_(status);}
	status = clWaitForEvents (1, &ev);		if (status != CL_SUCCESS)	{ cout << "\nclWaitForEventsh(1, &ev)="	 <<checkerror(status)  <<"\n"<<flush; exit_(status);}

											if(verbosity>0) {
												ss << (keyFrameCount*1000 + costVolCount);																							// Save buffers to file ###########
												DownloadAndSave_3Channel(imgmem, ss.str(), paths.at("imgmem"), width * height * sizeof(float)*3, baseImage_size,   CV_32FC3, 	false );
												DownloadAndSave(		 lomem,  ss.str(), paths.at("lomem"),  width * height * sizeof(float),   baseImage_size,   CV_32FC1, 	false , 8); // a little more than the num images in costvol.
												DownloadAndSave(		 himem,  ss.str(), paths.at("himem"),  width * height * sizeof(float),   baseImage_size,   CV_32FC1, 	false , 8); //params[LAYERS]
												DownloadAndSave(		 amem,   ss.str(), paths.at("amem"),   width * height * sizeof(float),   baseImage_size,   CV_32FC1, 	false , params[MAX_INV_DEPTH]);
												DownloadAndSave(		 dmem,   ss.str(), paths.at("dmem"),   width * height * sizeof(float),   baseImage_size,   CV_32FC1, 	false , params[MAX_INV_DEPTH]);
												if(verbosity>1) cout << "\ncostVolCount="<<costVolCount;
												cout << "\ncalcCostVol chk13_finished\n" << flush;
											}
}

void RunCL::cacheGValue2(cv::Mat &bgray, float theta)
{
																																			if(verbosity>0) cout<<"\ncacheGValue2_chk0"<<flush;
	params[THETA]	=  theta;
	cl_int status;
	cl_event writeEvt;
	status = clEnqueueWriteBuffer(m_queue,  param_buf, CL_FALSE, 0, 16 * sizeof(float), params, 0, NULL, &writeEvt);						// WriteBuffer param_buf ##########
	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: cacheGValue2_chk0\n" 			   <<endl; exit_(status);}
	status = clFlush(m_queue); 					if (status!=CL_SUCCESS){cout<<"\nclFlush status = " 				<< checkerror(status)<<"\n"<<flush;exit_(status);}
	status = waitForEventAndRelease(&writeEvt); if (status!=CL_SUCCESS){cout<<"\nwaitForEventAndRelease status = " 	<< checkerror(status)<<"\n"<<flush;exit_(status);}

	stringstream ss;
	ss << "cacheGValue2";
	cl_int 		res;
	cl_event 	ev1, ev2;
	size_t 		bgraySize_bytes = bgray.total() * bgray.elemSize();

	if (basegraymem == 0) {
		basegraymem = clCreateBuffer(m_context, CL_MEM_READ_ONLY , bgraySize_bytes, 0, &res);
		status 		= clEnqueueWriteBuffer(uload_queue, basegraymem, CL_FALSE, 0, width * height * sizeof(float), bgray.data, 0, NULL, &ev1);	// WriteBuffer basegraymem ##########
		if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; exit_(status);}
		clFlush(uload_queue); clFinish(uload_queue); if (res != CL_SUCCESS)	{ cout << "\nres = " << checkerror(res) <<"\n"<<flush; exit_(res);}
	}
	res = clSetKernelArg(cache3_kernel, 0, sizeof(cl_mem), &basegraymem);if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(cache3_kernel, 1, sizeof(cl_mem), &gxmem);		 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(cache3_kernel, 2, sizeof(cl_mem), &gymem);		 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(cache3_kernel, 3, sizeof(cl_mem), &g1mem);		 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(cache3_kernel, 4, sizeof(cl_mem), &param_buf);	 if(res!=CL_SUCCESS){cout<<"\nparam_buf res = "<<checkerror(res)<<"\n"<<flush;exit_(res);}

	clFlush(m_queue); clFinish(m_queue); if (res != CL_SUCCESS)	{ cout << "\nbefore cache3_kernel res = " << checkerror(res) <<"\n"<<flush; exit_(res);}
	res = clEnqueueNDRangeKernel(m_queue, cache3_kernel, 1, 0, &global_work_size, 0, 0, NULL, &ev2);										// run kernel CacheG3(..) ########
	clFlush(m_queue); clFinish(m_queue); if (res != CL_SUCCESS)	{ cout << "\nafter cache3_kerne res = " << checkerror(res) <<"\n"<<flush; exit_(res);}

	ss << (keyFrameCount);
											if(verbosity>1) {
												cout << "\n ss.str() = " << ss.str() << "\n" << flush;
												cout << "\n global_work_size2 = " << global_work_size << "\n" << flush;
												cout << "\nm_queue =  " << m_queue << flush;
												DownloadAndSave(basegraymem, ss.str(), paths.at("basegraymem"), width * height * sizeof(float) , baseImage_size, CV_32FC1, 		false /*true*/, 1);
												DownloadAndSave(gxmem,		 ss.str(), paths.at("gxmem"), 		width * height * sizeof(float) , baseImage_size, CV_32FC1, 		false /*true*/, 1);
												DownloadAndSave(gymem, 		 ss.str(), paths.at("gymem"), 		width * height * sizeof(float) , baseImage_size, CV_32FC1, 		false /*true*/, 1);
												DownloadAndSave(g1mem, 		 ss.str(), paths.at("g1mem"), 		width * height * sizeof(float) , baseImage_size, CV_32FC1, 		false /*true*/, 1);
											}
	keyFrameCount++;
																																			if(verbosity>0) cout<<"\ncacheGValue_finished"<<flush;
}

void RunCL::updateQD(float epsilon, float theta, float sigma_q, float sigma_d)
{
	params[EPSILON]			=  epsilon;
	params[SIGMA_Q]			=  sigma_q;
	params[SIGMA_D]			=  sigma_d;
	params[THETA]			=  theta;

	cl_int status, res;
	cl_event writeEvt, ev;
	status = clEnqueueWriteBuffer(uload_queue, param_buf, CL_FALSE, 0, 16 * sizeof(float), params, 0, NULL, &writeEvt); 					// WriteBuffer param_buf ##########
	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: updateQD_chk0\n"	<< endl; exit_(status);}
	status = clFlush(uload_queue); 						if (status != CL_SUCCESS)	{ cout << "\nclFlush status = " 				<< checkerror(status) <<"\n"<<flush; exit_(status);}
	status = waitForEventAndRelease(&writeEvt); 	if (status != CL_SUCCESS)	{ cout << "\nwaitForEventAndRelease status = " 	<< checkerror(status) <<"\n"<<flush; exit_(status);}

	stringstream ss;
	ss << "updateQD";
	res = clSetKernelArg(updateQD_kernel, 0, sizeof(cl_mem), &g1mem);	 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res) <<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateQD_kernel, 1, sizeof(cl_mem), &qmem);	 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res) <<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateQD_kernel, 2, sizeof(cl_mem), &dmem);	 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res) <<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateQD_kernel, 3, sizeof(cl_mem), &amem);	 if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res) <<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateQD_kernel, 4, sizeof(cl_mem), &param_buf);if(res!=CL_SUCCESS){cout<<"\nparam_buf res = "<<checkerror(res)<<"\n"<<flush;exit_(res);}

	res = clEnqueueNDRangeKernel(m_queue, updateQD_kernel, 1, 0, &global_work_size, 0, 0, NULL, &ev); 										// run updateQD_kernel  aka UpdateQD(..) ####
	if (res != CL_SUCCESS)	{ cout << "\nres = " << checkerror(res) <<"\n"<<flush; exit_(res);}
	status = clFlush(m_queue);				if (status!=CL_SUCCESS){cout<<"\nclFlush(m_queue) status ="  <<checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clWaitForEvents (1, &ev);		if (status!=CL_SUCCESS){cout<<"\nclWaitForEventsh(1, &ev) =" <<checkerror(status) <<"\n"<<flush; exit_(status);}

											if(verbosity>1){
												//size_t st  = width * height * sizeof(float);
												QDcount++;
												int this_count = count + QDcount;
												ss << this_count;
												cv::Size q_size( baseImage_size.width, 2* baseImage_size.height ); // 2x sized for qx and qy.
												DownloadAndSave(qmem,   ss.str(), paths.at("qmem"),    2*width * height * sizeof(float), q_size        , CV_32FC1, false , -1*params[MAX_INV_DEPTH]);
												DownloadAndSave(dmem,   ss.str(), paths.at("dmem"),    width * height * sizeof(float)  , baseImage_size, CV_32FC1, false , params[MAX_INV_DEPTH]);
												cout<<"\nRunCL::updateQD_chk3_finished\n"<<flush;
											}
}

void RunCL::updateA(float lambda, float theta)
{
	params[THETA]			=  theta;
	params[LAMBDA]			=  lambda;
	cl_int status;
	cl_event writeEvt;
	status = clEnqueueWriteBuffer(uload_queue,  param_buf, CL_FALSE, 0, 16 * sizeof(float), params, 0, NULL, &writeEvt);						// WriteBuffer param_buf ##########
	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; cout << "Error: \nRunCL::updateA_chk0\n" << endl; exit_(status);}
																																			else if(verbosity>0) {cout << "\nRunCL::updateA_chk0\t\tlayers="<< params[LAYERS] <<" \n" << flush;}
	status = clFlush(uload_queue); 					if (status != CL_SUCCESS)	{ cout << "\nclFlush status = " << status << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = waitForEventAndRelease(&writeEvt); if (status != CL_SUCCESS)	{ cout << "\nwaitForEventAndRelease status = "<<status<<checkerror(status)<<"\n"<<flush; exit_(status);}

	stringstream ss;
	ss << "updateA";
	cl_int res;
	cl_event ev;
	res = clSetKernelArg(updateA_kernel, 0, sizeof(cl_mem), &cdatabuf); if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateA_kernel, 1, sizeof(cl_mem), &amem);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateA_kernel, 2, sizeof(cl_mem), &dmem);		if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateA_kernel, 3, sizeof(cl_mem), &lomem);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateA_kernel, 4, sizeof(cl_mem), &himem);	if(res!=CL_SUCCESS){cout<<"\nres = "<<checkerror(res)<<"\n"<<flush;exit_(res);}
	res = clSetKernelArg(updateA_kernel, 5, sizeof(cl_mem), &param_buf);if(res!=CL_SUCCESS){cout<<"\nparam_buf res = "<<checkerror(res)<<"\n"<<flush;exit_(res);} // param_buf

	status = clFlush(m_queue); 				if (status != CL_SUCCESS)	{ cout << "\nclFlush(m_queue) status = " << checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clFinish(m_queue); 			if (status != CL_SUCCESS)	{ cout << "\nclFinish(m_queue)="<<status<<" "<<checkerror(status)<<"\n"<<flush; exit_(status);}

	res = clEnqueueNDRangeKernel(m_queue, updateA_kernel, 1, 0, &global_work_size, 0, 0, NULL, &ev); 										// run updateA_kernel  aka UpdateA(..) #####
	if (res != CL_SUCCESS)	{ cout << "\nres = " << checkerror(res) <<"\n"<<flush; exit_(res);}
	status = clFlush(m_queue);				if (status != CL_SUCCESS)	{ cout << "\nclFlush(m_queue) status = "<<status<<" "<< checkerror(status) <<"\n"<<flush; exit_(status);}
	status = clWaitForEvents (1, &ev);		if (status != CL_SUCCESS)	{ cout << "\nclWaitForEventsh(1, &ev)="	<<status<<" "<<checkerror(status)  <<"\n"<<flush; exit_(status);}

																																			if(verbosity>0) cout<<"\nRunCL::updateA_chk1"<<flush;
	count = keyFrameCount*1000000 + A_count*1000 + 999;
	A_count++;
											if(A_count%1==0 && verbosity>0){
												ss << count << "_theta"<<theta<<"_";
												cv::Size q_size( baseImage_size.width, 2* baseImage_size.height );
												DownloadAndSave(amem,   ss.str(), paths.at("amem"),    width * height * sizeof(float), baseImage_size, CV_32FC1,  false , params[MAX_INV_DEPTH]);
												DownloadAndSave(dmem,   ss.str(), paths.at("dmem"),    width * height * sizeof(float), baseImage_size, CV_32FC1,  false , params[MAX_INV_DEPTH]);
												DownloadAndSave(qmem,   ss.str(), paths.at("qmem"),  2*width * height * sizeof(float), q_size        , CV_32FC1,  false , 0.1 );
											}
																																			if(verbosity>0) cout<<"\nRunCL::updateA_chk2_finished,   params[MAX_INV_DEPTH]="<<params[MAX_INV_DEPTH]<<flush;
}

RunCL::~RunCL()
{
	cl_int status;
	status = clReleaseKernel(cost_kernel);      	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseKernel(cache3_kernel);		if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseKernel(updateQD_kernel);		if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseKernel(updateA_kernel);		if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }

	status = clReleaseProgram(m_program);			if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(m_queue);		if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(uload_queue);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(dload_queue);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseCommandQueue(track_queue);	if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
	status = clReleaseContext(m_context);			if (status != CL_SUCCESS)	{ cout << "\nstatus = " << checkerror(status) <<"\n"<<flush; }
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
