void
calc_approximate(const Mat1f &bt11_smooth, const Mat1f &bt11_clear, const Mat1b &land_mask, const Mat1b &invalid_mask, 
                 Mat1f &bt11_approx, int time_size)
{
	int x,y,z;
	double d1,d2,coeff;
	//placeholder
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
				for(z=0;z<time_size;z++){
					
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
					for(z=0;z<time_size;z++){
						bt11_approx(y,x,z)= coeff * bt11_smooth(y,x,z);
					}
				}
			}
		}
	}
}


void
approx_clear(const vector<string> &smooth_paths, const vector<string> &clear_paths, 
			 const vector<string> &original_paths,const Mat1b &land_mask, const Mat1b &invalid_mask, vector<string> &approx_paths, bool interp){	


    int time_size = smooth_paths.size();
    //int time_size = 100;
    //interpolate checks the n+1 index to ensure their derivative is small
	//approx is HEIGHT x WIDTH x time_size-1
    int approx_time = time_size-1;
    //int time_size = 10;
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG+SMOOTH_WINDOW_LAG;

    int i,j;
    int dims[3] = {HEIGHT,WIDTH,time_size};
    int dims_approx[3] = {HEIGHT,WIDTH,time_size-1};


    Mat1f approx(3,dims_approx);
    Mat1f clear_samples(3,dims);
    Mat1f smooth_samples(3,dims);
    Mat1b clear_masks(3,dims);

  	string approx_folder_loc = "data/approx";
  	string filename, save_loc;
    


    //generate save paths for approx values
    for(j=0;j<approx_time;j++){
    	filename = generate_filename(smooth_paths[j]);
    	save_loc = approx_folder_loc+filename;
    	approx_paths.push_back(save_loc);
    }


    printf("getting ice masks and cloud masks\n");
    for(j=0;j<time_size;j++){
    	read_mask(clear_paths[j+SMOOTH_WINDOW_LAG],clear_masks,j);
    	readgranule_oneband(original_paths[j+original_lag],clear_samples,j,"sea_surface_temperature");
		apply_mask_slice(clear_masks,clear_samples,j,true);
		readgranule_oneband(smooth_paths[j],smooth_samples,j,"sea_surface_temperature");
		printf("reading file %d for variable %s\n",j,"sea_surface_temperature");
    }

    printf("interpolating smooth samples\n");
    interpolate(smooth_samples, land_mask, invalid_mask, time_size, T_INTERP, PASS_THRESH, true, false);
    printf("interpolation completed\n");

    printf("calculating approximation for brightness_temperature_08um6\n");
 	calc_approximate(smooth_samples,clear_samples, land_mask,invalid_mask, approx,approx_time);
 	printf("finished approximation for brightness_temperature_08um6\n");
 	
 	//TODO WRITE THIS FUNCTION 
 	printf("Starting Collation on brightness_temperature_08um6\n");
 	collate_samples(clear_samples,approx, collated, Gamma,  collated_inds,land_mask,invalid_mask,interp);
    printf("Finished Collation on brightness_temperature_08um6\n");
    
    if(interp){
    	interpolate(collated,land_mask,invalid_mask,collated_interp_size,T_INTERP,false,false);
    	save_mat(collated_interp_paths, collated, "sea_surface_temperature",true);
    }
    else{
    	save_mat(collated_paths, collated, "sea_surface_temperature",true);
    }
    save_mat(approx_paths, approx, "sea_surface_temperature",true);

    
    
    /*
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
        save_approx(save_loc, save_slice,"brightness_temperature_10um4", false);
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
        save_approx(save_loc, save_slice,"brightness_temperature_11um2", false);
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
        save_approx(save_loc, save_slice,"brightness_temperature_12um3", false);
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
        save_approx(save_loc, save_slice,"sea_surface_temperature", false);
        approx_paths.push_back(save_loc.c_str());
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;
    */
}
