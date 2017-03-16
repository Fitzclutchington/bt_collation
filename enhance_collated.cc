//open all collated
//interpolate
//save new collated 
void
enhance_collated(const vector<string> &collated_paths, const vector<string> &approx_paths, const Mat1b &land_mask, const Mat1b &invalid_mask, bool interp)
{
	int y,x,t,i;
	int time_size = collated_paths.size();
	int mid = 3;
	int dims[3] = {HEIGHT, WIDTH,time_size};

	string folder_loc = "data/collated_interp";
	string save_loc, filename, temp_filename;

	vector<string> collated_files;
	vector<int> collated_inds;

	Mat1f bt08_collated(3,dims);
	Mat1f bt10_collated(3,dims);
	Mat1f bt11_collated(3,dims);
	Mat1f bt12_collated(3,dims);
	Mat1f sst_collated(3,dims);

	Mat1f bt08_slice(HEIGHT,WIDTH);
	Mat1f bt10_slice(HEIGHT,WIDTH);
	Mat1f bt11_slice(HEIGHT,WIDTH);
	Mat1f bt12_slice(HEIGHT,WIDTH);
	Mat1f sst_slice(HEIGHT,WIDTH);

	printf("reading collated files\n");
	for(t=0;t<collated_paths.size();t++){
		readgranule_fullbands(collated_paths[t].c_str(), bt11_collated,bt12_collated,bt08_collated,bt10_collated, sst_collated, t);
		printf("opened: %s\n",collated_paths[t].c_str());
	}
	printf("read all collated files\n");

	printf("starting interpolation of collated\n");
	interpolate(bt08_collated, land_mask, invalid_mask, time_size, T_COLL_INTERP, PASS_THRESH, false, false);
	interpolate(bt10_collated, land_mask, invalid_mask, time_size, T_COLL_INTERP, PASS_THRESH, false,false);
	interpolate(bt11_collated, land_mask, invalid_mask, time_size, T_COLL_INTERP, PASS_THRESH, false,false);
	interpolate(bt12_collated, land_mask, invalid_mask, time_size, T_COLL_INTERP, PASS_THRESH, false,false);
	interpolate(sst_collated, land_mask, invalid_mask, time_size, T_COLL_INTERP, PASS_THRESH, false,false);
	printf("finished interpolation of collated\n");

	if(interp){
		for(i=mid;i<approx_paths.size() - mid;i++){
	    	temp_filename = generate_filename(approx_paths[i]);
	    	if(temp_filename[11] == '0' && temp_filename[12] == '0'){
	    		collated_files.push_back(folder_loc+temp_filename);
	    		collated_inds.push_back(collated_files.size() - 1);
	    	}
	    	else if(collated_files.size()>0){
	    		collated_files.push_back(folder_loc+temp_filename);
	    	}

	    }

    for(i=0;i<collated_inds.size();i++){
    	printf("colalted at %d\n",collated_inds[i]);
    }
    /*
    for(i=mid;i<approx_paths.size() - mid;i++){
    	printf("file name = %s\n",collated_files[i].c_str());
    }
    */
    
	    printf("Starting interpolation of bt08\n");
	    int interp_dims[3] = {HEIGHT,WIDTH,collated_files.size()};
	    Mat1f collated_interp(3,interp_dims);

	    for(t=0;t<collated_inds.size();t++){
	    	for(y=0;y<HEIGHT;y++){
	    		for(x=0;x<WIDTH;x++){
	    			collated_interp(y,x,collated_inds[t]) = bt08_collated(y,x,t);
	    		}
	    	}
	    }

	    interpolate(collated_interp, land_mask, invalid_mask, time_size, PASS_THRESH, INTERP_DIST, false, false);
	    save_mat(collated_files, collated_interp, "brightness_temperature_08um6",true);
	    printf("finished interpolation of bt08\n");

	    for(t=0;t<collated_inds.size();t++){
	    	for(y=0;y<HEIGHT;y++){
	    		for(x=0;x<WIDTH;x++){
	    			collated_interp(y,x,collated_inds[t]) = bt10_collated(y,x,t);
	    		}
	    	}
	    }

	    interpolate(collated_interp, land_mask, invalid_mask, time_size, PASS_THRESH, INTERP_DIST, false, false);
	    save_mat(collated_files, collated_interp, "brightness_temperature_10um4",false);

	    for(t=0;t<collated_inds.size();t++){
	    	for(y=0;y<HEIGHT;y++){
	    		for(x=0;x<WIDTH;x++){
	    			collated_interp(y,x,collated_inds[t]) = bt11_collated(y,x,t);
	    		}
	    	}
	    }

	    interpolate(collated_interp, land_mask, invalid_mask, time_size, PASS_THRESH, INTERP_DIST, false, false);
	    save_mat(collated_files, collated_interp, "brightness_temperature_11um2",false);

	     for(t=0;t<collated_inds.size();t++){
	    	for(y=0;y<HEIGHT;y++){
	    		for(x=0;x<WIDTH;x++){
	    			collated_interp(y,x,collated_inds[t]) = bt12_collated(y,x,t);
	    		}
	    	}
	    }

	    interpolate(collated_interp, land_mask, invalid_mask, time_size, PASS_THRESH, INTERP_DIST, false, false);
	    save_mat(collated_files, collated_interp, "brightness_temperature_12um3",false);

	     for(t=0;t<collated_inds.size();t++){
	    	for(y=0;y<HEIGHT;y++){
	    		for(x=0;x<WIDTH;x++){
	    			collated_interp(y,x,collated_inds[t]) = sst_collated(y,x,t);
	    		}
	    	}
	    }

	    interpolate(collated_interp, land_mask, invalid_mask, time_size, PASS_THRESH, INTERP_DIST, false, false);
	    save_mat(collated_files, collated_interp, "sea_surface_temperature",false);
	}
	else{
		for(i=0;i<collated_paths.size();i++){
			filename = generate_filename(collated_paths[i]);
			collated_files.push_back(folder_loc+filename);
			printf("collated file = %s\n",collated_files[i].c_str());
		}

		save_mat(collated_files,bt08_collated,"brightness_temperature_08um6",true);
		save_mat(collated_files,bt10_collated,"brightness_temperature_10um4",false);
		save_mat(collated_files,bt11_collated,"brightness_temperature_11um2",false);
		save_mat(collated_files,bt12_collated,"brightness_temperature_12um3",false);
		save_mat(collated_files,sst_collated,"sea_surface_temperature",false);
	}
    /*
	for(t=0;t<collated_paths.size();t++){

		for(y=0;y<HEIGHT;y++){
			for(x=0;x<WIDTH;x++){
				bt08_slice(y,x) = bt08_collated(y,x,t);
				bt10_slice(y,x) = bt10_collated(y,x,t);
				bt11_slice(y,x) = bt11_collated(y,x,t);
				bt12_slice(y,x) = bt12_collated(y,x,t);
				sst_slice(y,x) = sst_collated(y,x,t);
			}
		}
		

		filename = generate_filename(collated_paths[t]);
		save_loc = folder_loc + filename;
		save_test_nc_fullbands(bt08_slice, bt10_slice, bt11_slice, bt12_slice, sst_slice, save_loc.c_str());
	}
	*/
	
}	