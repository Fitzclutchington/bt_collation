void
calc_approximate(const Mat1f &bt11_smooth, const Mat1f &bt11_clear, const Mat1b &land_mask, const Mat1b &invalid_mask, 
                 Mat1f &bt11_approx, int * dims)
{
	int x,y,z,i;
	double d1,d2,coeff;
	//placeholder
	int time = dims[2];
    int count=0;
	vector<float> data_cl2;
	vector<float> smooth;
	vector<int> clear_inds;
	//Mat1f diffs(HEIGHT,WIDTH);
	bt11_approx.setTo(NAN);
	for(y=0;y<HEIGHT;y++){
		for(x=0;x<WIDTH;x++){
			if(land_mask(y,x) == 0 && invalid_mask(y,x) == 0 ){
                count = 0;
				d1=0;d2 =0;
				data_cl2.clear();
				smooth.clear();
				clear_inds.clear();
				for(z=0;z<time;z++){
					
					if(!std::isnan(bt11_smooth(y,x,z)) && !std::isnan(bt11_clear(y,x,z))){
					    //if(!std::isnan(bt11_clear(y,x,z)) || (fabs(bt11_smooth(y,x,z) - sst(y,x,z))< T_SST_SMOOTH)){

						d1 += bt11_clear(y,x,z)*bt11_smooth(y,x,z);
						d2 += bt11_smooth(y,x,z)*bt11_smooth(y,x,z);
                        count++;
                    			
					}					
				}
			
				//diffs(y,x) = d2-d1;
				if(d2 != 0 && count > MIN_POINTS){
					coeff = d1/d2;
					for(z=0;z<time;z++){
						bt11_approx(y,x,z)= coeff * bt11_smooth(y,x,z);
					}
				}
			}
		}
	}
}



void
generate_new_clear(const Mat1f &diffs,const Mat1f &original_samples,Mat1f &clear_samples,int &time_size)
{
	int x,y,t;
	float delta = 0.3;
	for(y=0;y<HEIGHT;y++){
		for(x=0;x<WIDTH;x++){
			for(t=0;t<time_size;t++){
				clear_samples(y,x,t) = NAN;
				if(diffs(y,x,t) < delta && !std::isnan(diffs(y,x,t))){
					clear_samples(y,x,t) = original_samples(y,x,t);
				}
			}
		}
	}
}

void
approximate(const vector<string> &original_paths,const vector<string> &smooth_paths,const Mat1b &clear_masks, const Mat1b &land_mask, 
	        const Mat1b &invalid_mask, Mat1f &clear_samples, Mat1f &smooth_samples, const string variable,Mat1f &approx,int time_size)
{
	int j;
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG+SMOOTH_WINDOW_LAG;
	int dims[3] = {HEIGHT,WIDTH,time_size-1};


    //TODO: Make function so only interp and clear are returned, don't need smooth
	for(j=0;j<time_size;j++){

		readgranule_oneband(original_paths[j+original_lag],clear_samples,j,variable);
		apply_mask_slice(clear_masks,clear_samples,j,true);
		readgranule_oneband(smooth_paths[j],smooth_samples,j,variable);
		printf("reading file %d for variable %s\n",j,variable.c_str());
	}

	printf("starting approximation\n");

    printf("interpolating smooth samples\n");
    interpolate(smooth_samples, land_mask, invalid_mask, time_size, T_INTERP, PASS_THRESH, true, false);
    printf("interpolation completed\n");

 	calc_approximate(smooth_samples,clear_samples, land_mask,invalid_mask, approx,dims); 

    printf("finished approximation\n");
}

void
approx_clear(const vector<string> &smooth_paths, const vector<string> &clear_paths, 
			 const vector<string> &original_paths,const Mat1b &land_mask, const Mat1b &invalid_mask, vector<string> &approx_paths){	
	// Compute approximation

    int y,x;
    int time_size = smooth_paths.size();
    //int time_size = 10;
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG+SMOOTH_WINDOW_LAG;

    int i,j;
    int file_count =0;
    int start_ind = 0;
    int dims[3] = {HEIGHT,WIDTH,time_size};
    int dims_approx[3] = {HEIGHT,WIDTH,time_size-1};
    int inds[time_size];

    Mat1f save_slice(HEIGHT,WIDTH);

    Mat1f approx(3,dims_approx);
    Mat1f clear_samples(3,dims);
    Mat1f smooth_samples(3,dims);
    Mat1b clear_masks(3,dims);
    //Mat1f diffs(3,dims);

  	compute_indicies(inds, start_ind, time_size);
  	string folder_loc = "data/approx";
  	string filename, save_loc;
    //smooth sample files
    
    //TODO combine these in a function to save space
    printf("getting ice masks and cloud masks\n");
    for(j=0;j<time_size;j++){
    	read_mask(clear_paths[j+SMOOTH_WINDOW_LAG],clear_masks,j);
    	inds[j] = j;
    }
    printf("generated ice and cloud mask\n");

    //TODO: append to file in functions to save space, don't need all bands saved at one time
  	approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples,smooth_samples, "brightness_temperature_08um6",approx, time_size);

    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
    	printf("filename = %s\n",filename);
    	file_count++;

    	save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_and_update(save_loc, save_slice,"brightness_temperature_08um6", true);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count=0;
    
    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "brightness_temperature_10um4",approx, time_size);

    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_and_update(save_loc, save_slice,"brightness_temperature_10um4", false);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;

    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "brightness_temperature_11um2",approx, time_size);
    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_and_update(save_loc, save_slice,"brightness_temperature_11um2", false);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;

    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "brightness_temperature_12um3",approx, time_size);
    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_and_update(save_loc, save_slice,"brightness_temperature_12um3", false);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;
    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "sea_surface_temperature",approx, time_size);
    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_and_update(save_loc, save_slice,"sea_surface_temperature", false);
        approx_paths.push_back(save_loc.c_str());
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;
}

