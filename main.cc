#include <netcdf.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdarg.h>
#include <cmath>
#include <iostream>
#include <fstream> 
#include <Eigen/Dense>

using namespace std;
using namespace cv;
using namespace Eigen;
#include "parameters.cc"
#include "io.cc"
#include "calc.cc"
#include "mask.cc"
#include "filter.cc"
#include "smoothing.cc"
#include "interpolate.cc"
#include "approximate.cc"
#include "subsample.cc"
#include "second_pass.cc"
#include "enhance_collated.cc"
//#include "generateLUT.cc"
#define NDEBUG


int
main(int argc, char *argv[])
{
    int i;
    if(argc < 2){
        eprintf("Usage: ./ahil2c <granule_list>");
    }

    Mat1b land_mask(HEIGHT,WIDTH);
    Mat1b invalid_mask(HEIGHT,WIDTH);
    Mat1b border_mask(HEIGHT,WIDTH);

    vector<string> original_paths;
    vector<string> second_pass_paths;
    vector<string> clear_paths;
    vector<string> smooth_paths;
    vector<string> approx_paths;
    vector<string> collated_paths;

    get_file_names(argv[1],original_paths);
    //string reference_file = argv[2];

    //function to get land mask and invalid(corners) mask
    get_l2pmask(original_paths[0].c_str(),land_mask,invalid_mask);
    //get_lats(original_paths[0], lats);
    get_landborders(land_mask, border_mask, LAND_KERNEL);
    //SAVENC(border_mask);
    printf("read land mask\n");
    string filename;
    string clearpath, smoothpath, approxpath;
        
    
    //for debugging purposes
    for(i =0;i<222;i++){
        filename = generate_filename(original_paths[i+FILTER_WINDOW_LAG]);
        clearpath = "data/clear" + filename;
        clear_paths.push_back(clearpath.c_str());
    }
    
    /*
    for(i=0;i<198;i++){
        filename = generate_filename(original_paths[i+FILTER_WINDOW_LAG+SECOND_PASS_LAG]);
        clearpath = "data/pass2" + filename;
        second_pass_paths.push_back(clearpath.c_str());
    }
    */
    /*
    setupLUT(second_pass_paths, original_paths, reference_file);
    */

    /*
    printf("starting cloud filter\n");
    filter_clouds(original_paths, land_mask,invalid_mask, border_mask, clear_paths);
    printf("Cloud filter completed for %d granules\n",clear_paths.size());  
    */
    
    
    second_pass(clear_paths, original_paths, second_pass_paths, land_mask, invalid_mask);
    
    /*
    for(i =0;i<180;i++){
        filename = generate_filename(original_paths[i+FILTER_WINDOW_LAG+SMOOTH_WINDOW_LAG+SECOND_PASS_LAG]);
        smoothpath = "data/smooth" + filename;
        smooth_paths.push_back(smoothpath.c_str());
    }
    */
      

    
    printf("starting smoothing operation\n");
    smooth_samples(second_pass_paths, original_paths, land_mask, invalid_mask, smooth_paths);
    printf("Smoothing completed for %d granules\n",smooth_paths.size());
    
   
    
    
    /*
    for(i =0;i<179;i++){
        filename = generate_filename(smooth_paths[i]);
        approxpath = "data/approx" + filename;
        approx_paths.push_back(approxpath.c_str());
    }    
    */
    
    
    printf("starting approximation\n");
    approx_clear(smooth_paths, second_pass_paths, original_paths,land_mask, invalid_mask,approx_paths);
    printf("SApproximation completed for %d granules\n",approx_paths.size());
      
    
    
    printf("starting subsample\n");
    subsample(approx_paths,original_paths, smooth_paths, collated_paths, second_pass_paths, land_mask,invalid_mask);
    printf("Subsample completed for %d granules\n",approx_paths.size());
    
    
    /*
    for(i =3;i<176;i++){
        filename = generate_filename(smooth_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            approxpath = "data/collated" + filename;
            collated_paths.push_back(approxpath.c_str());
        }
    }  
    */

    printf("starting subsample\n");
    enhance_collated(collated_paths, land_mask, invalid_mask);
    printf("Subsample completed for %d granules\n",approx_paths.size());
    
}
