//open all collated
//interpolate
//save new collated 
void
enhance_collated(const vector<string> &collated_paths, const Mat1b &land_mask, const Mat1b &invalid_mask)
{
	int y,x,t;
	int time_size = collated_paths.size();
	int dims[3] = {HEIGHT, WIDTH,time_size};

	string folder_loc = "data/collated_interp";
	string save_loc, filename;

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
	interpolate(bt08_collated, land_mask, invalid_mask, time_size, T_INTERP, false, false);
	interpolate(bt10_collated, land_mask, invalid_mask, time_size, T_INTERP, false,false);
	interpolate(bt11_collated, land_mask, invalid_mask, time_size, T_INTERP, false,false);
	interpolate(bt12_collated, land_mask, invalid_mask, time_size, T_INTERP, false,false);
	interpolate(sst_collated, land_mask, invalid_mask, time_size, T_INTERP, false,false);
	printf("finished interpolation of collated\n");

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
		/*
		update_variable(collated_paths[t], bt08_slice,"brightness_temperature_08um6");
		update_variable(collated_paths[t], bt10_slice,"brightness_temperature_10um4");
		update_variable(collated_paths[t], bt11_slice,"brightness_temperature_11um2");
		update_variable(collated_paths[t], bt12_slice,"brightness_temperature_12um3");
		update_variable(collated_paths[t], sst_slice,"sea_surface_temperature");
		printf("generated path for collated %d\n",t);
		*/
		filename = generate_filename(collated_paths[t]);
		save_loc = folder_loc + filename;
		save_test_nc_fullbands(bt08_slice, bt10_slice, bt11_slice, bt12_slice, sst_slice, save_loc.c_str());
	}
}	