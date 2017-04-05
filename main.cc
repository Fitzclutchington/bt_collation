#include <netcdf.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdarg.h>
#include <cmath>
#include <iostream>
#include <fstream> 
#include <Eigen/Dense>
#include <chrono>

using namespace std;
using namespace cv;
using namespace Eigen;
#include "hermite.h"
#include "parameters.cc"
#include "io.cc"
#include "calc.cc"
#include "mask.cc"
#include "filter.cc"
#include "smooth_no_padding.cc"
#include "interpolate.cc"
#include "second_pass.cc"
#include "collated_fit_no_padding.cc"
#define NDEBUG


int
main(int argc, char *argv[])
{


    auto start = std::chrono::system_clock::now();
    int window[3] = {CLEAR_SPATIAL_SMOOTH, CLEAR_SPATIAL_SMOOTH, SMOOTH_WINDOW_LAG};

    if(argc < 4){
        eprintf("Usage: ./ahil2c <granule_list> <reference_file> <interp>");
    }

    bool interp = true;
    if(strcmp(argv[3], "--interpolate") != 0){
        interp = false;
    }
    else{
        interp = true;
        printf("interpolating\n");
    }

    Mat1b land_mask(HEIGHT,WIDTH);
    Mat1b invalid_mask(HEIGHT,WIDTH);
    Mat1b border_mask(HEIGHT,WIDTH);
    Mat1b l2p_mask(HEIGHT,WIDTH);
    Mat1b ice_mask(HEIGHT,WIDTH);

    vector<string> original_paths;
    vector<string> second_pass_paths;
    vector<string> clear_paths;
    vector<string> smooth_paths;
    vector<string> approx_paths;
    vector<string> collated_paths;

    get_file_names(argv[1],original_paths);
    string reference_file = argv[2];

    //function to get land mask and invalid(corners) mask
    get_l2pmask(original_paths[0].c_str(),land_mask,invalid_mask);
    
    //get ice_mask
    get_icemask(reference_file, ice_mask);


    //get_lats(original_paths[0], lats);
    get_landborders(land_mask, border_mask, LAND_KERNEL);
    
    combine_l2pmasks(land_mask,invalid_mask, ice_mask, l2p_mask);
    printf("read land mask\n");
    land_mask.release();
    invalid_mask.release();
    ice_mask.release();
    //SAVENC(border_mask);
    
    //SAVENC(l2p_mask);
    string filename;
    string clearpath, smoothpath, approxpath;
        
    /*
    //for debugging purposes
    int i;
    for(i =0;i<203;i++){
        filename = generate_filename(original_paths[i+FILTER_WINDOW_LAG]);
        clearpath = "data/clear" + filename;
        clear_paths.push_back(clearpath.c_str());
    }
    
    
    for(i=0;i<179;i++){
        filename = generate_filename(original_paths[i+FILTER_WINDOW_LAG+SECOND_PASS_LAG]);
        clearpath = "data/pass2" + filename;
        second_pass_paths.push_back(clearpath.c_str());
    }

    for(i =0;i<179;i++){
        filename = generate_filename(original_paths[i+FILTER_WINDOW_LAG+SECOND_PASS_LAG]);
        smoothpath = "data/smooth" + filename;
        smooth_paths.push_back(smoothpath.c_str());
    }
    

    */
    
    
    printf("starting cloud filter\n");


    filter_clouds(original_paths, l2p_mask, border_mask, clear_paths, reference_file);
    

    int clear_size = clear_paths.size();
    printf("Cloud filter completed for %d granules\n",clear_size);      
    
    
    printf("starting second pass\n");

    second_pass(clear_paths, original_paths, second_pass_paths, l2p_mask);
    

    printf("finished second pass\n");
    
    
    printf("starting smoothing operation\n");
    smooth_samples(second_pass_paths, original_paths, l2p_mask, smooth_paths, reference_file, window, false);
    int smooth_size = smooth_paths.size();
    printf("Smoothing completed for %d granules\n",smooth_size);
    
    
    
    //approx_clear(smooth_paths,second_pass_paths, original_paths,land_mask, invalid_mask, approx_paths,reference_file,interp);
    printf("starting approximation and collation\n");

    approx_clear(smooth_paths, second_pass_paths, original_paths,l2p_mask, reference_file, interp);
    printf("finished approximation and collation\n");

    auto end = std::chrono::system_clock::now();
    auto elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to compute code = " << elapsed.count() << " seconds\n";
    
}
