
void
smoothing3D(const Mat1f &masked_data, const Mat1f &reference, const Mat1b &l2p_mask, 
            Mat1f &smooth_output, int *window, int cur_ind)
{
    int x,y,t,i,j,k;

    int w = 50;



    int t_len = 2*SMOOTH_WINDOW_LAG+1;
    int w_len = 2*w+1;
    int threshold;
    double left_sum, sum,  val, avg,  count, left_count;

    Mat1f time_sum(HEIGHT,WIDTH);
    Mat1f time_count(HEIGHT,WIDTH);

    threshold = PR_CLEAR * (2*window[0]+1)*(2*window[1]+1)*(2*window[2]+1);
    //threshold =1;
    //printf("%d \n",threshold);
    t = SMOOTH_WINDOW_LAG;
    printf("starting loop\n");

    //windowed_nanmean_3d(bt08, land_mask, invalid_mask, bt08_3dsmooth, window, threshold);
    //windowed_nanmean_3d(bt10, land_mask, invalid_mask, bt10_3dsmooth, window, threshold);
    //windowed_nanmean_3d(bt11, land_mask, invalid_mask, bt11_3dsmooth, window, threshold);
    //windowed_nanmean_3d(bt12, land_mask, invalid_mask, bt12_3dsmooth, window, threshold);
    //subtract ref from sst ( sst - ref)
    
    windowed_nanmean_3d(sst, land_mask, invalid_mask,  sst_3dsmooth, window, threshold);
    //add back reference to sst_3dsmooth ( sst_3dsmooth + ref)
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            sst_3dsmooth(y,x) = sst_3dsmooth(y,x) + reference(y,x);
        }
    }
    
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// smooth_samples:
//
// Input:
// vector<string> clear_paths - folder locations of cloud mask files 
// vector<string> original_paths- folder locations of original satellite data files
// Mat1b l2p_mask - a HEIGHT x WIDTH matrix where 2555 is a land or invalid pixel and 0 is a valid pixel 
// vector<string> smooth_paths - vector where smooth save location file names will be placed 
// string ref_file - location of reference file
// int *window - array containing the size of the full smoothing window lags {x_lag,y_lag,time_lag}
// bool full - if true compute on all bands, else compute only on sst
//
// This function is a controller for the smoothing operation.  It opens the original satellite data files
// and applies the computed mask from the cloud filtering step. 
//
void
smooth_samples(const vector<string> clear_paths, const vector<string> original_paths,const Mat1b &l2p_mask, 
               vector<string> &smooth_paths, const string ref_file, int *window, bool full){
    
    int i,y,x,t;
    int sample_size = clear_paths.size(); // number of clear files
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG; // lag to match original files to mask
    int dims[3] = {HEIGHT,WIDTH,sample_size}; // dimensions of array containing masked data
    string folder_loc = "data/smooth"; // folder to save output files
    string filename, save_loc;


    Mat1b clear_mask(HEIGHT,WIDTH); 
    
    Mat1f masked_data(3,dims);
    Mat1f smooth_output(HEIGHT, WIDTH); // matrix of smooth vals to save to file
    Mat1f reference(HEIGHT,WIDTH);

    // get the reference data
    get_var(ref_file,reference,"sst_reynolds");

    // always compute on sst
    // open sst and apply mask
    for(i = 0; i < sample_size; ++i){

        //open mask file
        read_mask(clear_paths[i],clear_mask,-1);
        //open original data file and place into buffer
        readgranule_oneband(original_paths[i+original_lag],masked_data,i,"sea_surface_temperature");
        //apply mask to current buffer slice
        apply_mask_slice(clear_mask,masked_data,i,false);

        //subtract reference from granule
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                sst(y,x,i) = sst(y,x,i) - reference(y,x);
            }
        }
        printf("read file %s\n",original_paths[i+original_lag].c_str());
    }

    // smooth masked_data matrix
    // run smoothing function given index
    // if index does not fill window, ensure window dimensions are changed
    // change window to {y_lag,x_lag,time_front,time_back}
    // user iterative process (front sum, front count, back sum, back count)
    for(i = 0; i < sample_size; ++i){
        //determine what window needs to be used

    }
        
}