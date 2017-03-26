
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
    
    int i,j,y,x,t,cur_ind;
    int threshold;
    float val;
    int sample_size = clear_paths.size(); // number of clear files
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG; // lag to match original files to mask
    int time_lag = window[2];
    int dims[3] = {HEIGHT,WIDTH,sample_size}; // dimensions of array containing masked data
    string folder_loc = "data/smooth"; // folder to save output files
    string filename, save_loc;


    Mat1b clear_mask(HEIGHT,WIDTH); 
    
    Mat1f masked_data(3,dims);
    Mat1f smooth_output(HEIGHT, WIDTH); // matrix of smooth vals to save to file
    Mat1f reference(HEIGHT,WIDTH);

    Mat1f time_sum(HEIGHT,WIDTH);
    time_sum.setTo(0);
    
    Mat1f time_count(HEIGHT,WIDTH);
    time_count.setTo(0);
    
    Mat1f left_time_sum(HEIGHT,WIDTH);
    Mat1f left_time_count(HEIGHT,WIDTH);
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
                masked_data(y,x,i) = masked_data(y,x,i) - reference(y,x);
            }
        }

        // save locations for smooth files
        filename = generate_filename(clear_paths[i]);
        smooth_paths.push_back(folder_loc+filename);
        //printf("%s\n",smooth_paths[i].c_str());
        printf("read file %s\n",original_paths[i+original_lag].c_str());
    }


    //determine first time sum and time count
    for(t = 0; t < SMOOTH_WINDOW_LAG + 1; ++t){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                val = masked_data(y,x,t);
                if(std::isfinite(val)){
                    time_count(y,x)++;
                    time_sum(y,x) += val;
                }
            }
        }
    }

    threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * time_lag;

    // compute first smoothb
    windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold); // rewrite this function to only work with time sum
    // add back reference
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            smooth_output(y,x) = smooth_output(y,x) + reference(y,x);
        }
    }
    save_and_update(smooth_paths[0],smooth_output,"sea_surface_temperature",true);


    //compute smooth before there arent enough granules to fill entire window
    for(cur_ind = 1; cur_ind < time_lag+1; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "right");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag+cur_ind);
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                smooth_output(y,x) = smooth_output(y,x) + reference(y,x);
            }
        }
        save_and_update(smooth_paths[cur_ind],smooth_output,"sea_surface_temperature",true); 
        printf("smoothed file %s\n",smooth_paths[cur_ind].c_str());
    }   

    //compute in interval where entire window will be full
    for(; cur_ind < sample_size - time_lag; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "both");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (2*time_lag +1);
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                smooth_output(y,x) = smooth_output(y,x) + reference(y,x);
            }
        }
        save_and_update(smooth_paths[cur_ind],smooth_output,"sea_surface_temperature",true); 
        printf("smoothed file %s\n",smooth_paths[cur_ind].c_str());
    }

    //compute end interval
    for(; cur_ind < sample_size; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "left");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag + sample_size - cur_ind); // i counts how many granules to ignore as window shrinks
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                smooth_output(y,x) = smooth_output(y,x) + reference(y,x);
            }
        }
        save_and_update(smooth_paths[cur_ind],smooth_output,"sea_surface_temperature",true); 
        printf("smoothed file %s\n",smooth_paths[cur_ind].c_str());
        
    }

    if(full){
        vector<string> band_names;
        band_names.push_back("brightness_temperature_08um6");
        band_names.push_back("brightness_temperature_10um4");
        band_names.push_back("brightness_temperature_11um2");
        band_names.push_back("brightness_temperature_12um3");
        int total_bands = band_names.size();
        
        for(j = 0; j < total_bands; ++j){
            time_sum.setTo(0);
            time_count.setTo(0);

            for(i = 0; i < sample_size; ++i){

                //open mask file
                read_mask(clear_paths[i],clear_mask,-1);
                //open original data file and place into buffer
                readgranule_oneband(original_paths[i+original_lag],masked_data,i,band_names[j]);
                //apply mask to current buffer slice
                apply_mask_slice(clear_mask,masked_data,i,false);

                //printf("%s\n",smooth_paths[i].c_str());
                printf("read file %s\n",original_paths[i+original_lag].c_str());
            }


            //determine first time sum and time count
            for(t = 0; t < SMOOTH_WINDOW_LAG + 1; ++t){
                for(y = 0; y < HEIGHT; ++y){
                    for(x = 0; x < WIDTH; ++x){
                        val = masked_data(y,x,t);
                        if(std::isfinite(val)){
                            time_count(y,x)++;
                            time_sum(y,x) += val;
                        }
                    }
                }
            }

            threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * time_lag;

            // compute first smoothb
            windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold); // rewrite this function to only work with time sum

            save_and_update(smooth_paths[0],smooth_output,band_names[j],false);


            //compute smooth before there arent enough granules to fill entire window
            for(cur_ind = 1; cur_ind < time_lag+1; ++cur_ind){
                update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "right");

                threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag+cur_ind);
                windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  

                save_and_update(smooth_paths[cur_ind],smooth_output,band_names[j],false); 
                printf("smoothed file %s\n",smooth_paths[cur_ind].c_str());
            }   

            //compute in interval where entire window will be full
            for(; cur_ind < sample_size - time_lag; ++cur_ind){
                update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "both");

                threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (2*time_lag +1);
                windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
       
                save_and_update(smooth_paths[cur_ind],smooth_output,band_names[j],false); 
                printf("smoothed file %s\n",smooth_paths[cur_ind].c_str());
            }

            //compute end interval
            for(; cur_ind < sample_size; ++cur_ind){
                update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "left");

                threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag + sample_size - cur_ind); // i counts how many granules to ignore as window shrinks
                windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  

                save_and_update(smooth_paths[cur_ind],smooth_output,band_names[j],false); 
                printf("smoothed file %s\n",smooth_paths[cur_ind].c_str());
                
            }       
        }    
    }
}

void
smooth_samples_collated(const Mat1f &collated_interp, Mat1f &collated_smooth,const Mat1f &l2p_mask,
                        const string ref_file,const int* window,const int sample_size)
{
    int i,y,x,t,cur_ind;
    int threshold;
    float val;
    
    int time_lag = window[2];
    int dims[3] = {HEIGHT,WIDTH,sample_size}; // dimensions of array containing masked data


    
    Mat1f smooth_output(HEIGHT, WIDTH); // matrix of smooth vals to save to file
    Mat1f reference(HEIGHT,WIDTH);

    Mat1f masked_data(3,dims);

    Mat1f time_sum(HEIGHT,WIDTH);
    time_sum.setTo(0);
    
    Mat1f time_count(HEIGHT,WIDTH);
    time_count.setTo(0);
    
    Mat1f left_time_sum(HEIGHT,WIDTH);
    Mat1f left_time_count(HEIGHT,WIDTH);

    // get the reference data
    get_var(ref_file,reference,"sst_reynolds");

    // always compute on sst
    // open sst and apply mask
    for(i = 0; i < sample_size; ++i){

        //subtract reference from granule
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                masked_data(y,x,i) = collated_interp(y,x,i) - reference(y,x);
            }
        }
    }

    printf("subtracted reference\n");

    //determine first time sum and time count
    for(t = 0; t < time_lag + 1; ++t){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                val = masked_data(y,x,t);
                if(std::isfinite(val)){
                    time_count(y,x)++;
                    time_sum(y,x) += val;
                }
            }
        }
    }
    printf("generated time sum\n");


    threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * time_lag;

    // compute first smoothb
    windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold); // rewrite this function to only work with time sum
    // add back reference
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            collated_smooth(y,x,0) = smooth_output(y,x) + reference(y,x);
        }
    }

    printf("generated first colalted smooth\n");
    //compute smooth before there arent enough granules to fill entire window
    for(cur_ind = 1; cur_ind < time_lag+1; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "right");
        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag+cur_ind);
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_smooth(y,x,cur_ind) = smooth_output(y,x) + reference(y,x);
            }
        } 
    }   

    printf("generated first colalted smooth interval\n");

    //compute in interval where entire window will be full
    for(; cur_ind < sample_size - time_lag; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "both");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (2*time_lag +1);
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_smooth(y,x,cur_ind) = smooth_output(y,x) + reference(y,x);
            }
        }
    }
    printf("generated middle colalted smooth\n");
    //compute end interval
    for(; cur_ind < sample_size; ++cur_ind){
        update_sums(time_count,time_sum,masked_data, cur_ind, time_lag, "left");

        threshold = PR_CLEAR *(2*window[0] + 1) * (2*window[1] + 1) * (time_lag + sample_size - cur_ind); // i counts how many granules to ignore as window shrinks
        windowed_nanmean_3d(time_count, time_sum, l2p_mask, smooth_output, window, threshold);  
        // add back reference
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_smooth(y,x,cur_ind) = smooth_output(y,x) + reference(y,x);
            }
        }         
    }
    printf("generated final colalted smooth\n");

    masked_data.release();
}