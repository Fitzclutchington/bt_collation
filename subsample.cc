
void
calc_subsample(const Mat1f &approx_samples, const Mat1f &approx_interp, const Mat1f &clear_samples, const Mat1f &samples, const Mat1b &land_mask, const Mat1b &invalid_mask, 
		       Mat1f &final_sample, Mat1b &clear, Mat1b &approx, int time_size, int ind, int mid)
{	
	float sum, collated_val, max_val, approx_val;
	int count_clear;
	int count_approx;
	int y,x,t;

	for(y=0;y<HEIGHT;y++){
			for(x=0;x<WIDTH;x++){
				final_sample(y,x) = NAN;
				clear(y,x) = 0;
				approx(y,x) = 0;
				count_approx =0;
				count_clear=0;
				sum = 0;
				if(land_mask(y,x) == 0 && invalid_mask(y,x) == 0){
					
					max_val = samples(y,x,0);
					for(t=0;t<time_size;t++){
						if(max_val < samples(y,x,t)){
							max_val = samples(y,x,t);
						}
						
						if(std::isnan(clear_samples(y,x,t))){
							if(!std::isnan(approx_samples(y,x,t + ind - mid))){
								sum+= MU_APPROX*approx_samples(y,x,t + ind - mid);
								count_approx+=1;
							}
						}
						else{
							sum+=MU_CLEAR*clear_samples(y,x,t);
							count_clear+=1;
						}
						
					}
					
					if(count_clear > 0 || count_approx > 0){
						collated_val = sum/(MU_CLEAR*count_clear +MU_APPROX*count_approx);
						approx_val = approx_interp(y,x,ind);
						if((max_val - collated_val < T_COLLATED || max_val - collated_val < 0) && (approx_val-collated_val < T_COLLATED || approx_val-collated_val < 0 )){
							final_sample(y,x) = collated_val;
							clear(y,x) = count_clear;
							approx(y,x) = count_approx;
						}
					}
				}
			}
		}
}


void
calc_subsample_3d(const Mat1f &approx_samples, const Mat1f &approx_interp, const Mat1f &clear_samples, const Mat1f &samples, const Mat1b &land_mask, 
				  const Mat1b &invalid_mask, Mat1f &final_sample, Mat1b &clear, Mat1b &approx, int time_size, int ind, int mid)
{	
	float sum, collated_val, max_val;
	int count_clear,count_approx;
	int y,x,t,i,j;
	int w = 1;
    
	final_sample.setTo(NAN);
    clear.setTo(0);
    approx.setTo(0);

	for(y=w;y<HEIGHT-w;y++){
			for(x=w;x<WIDTH-w;x++){

				count_approx =0;
				count_clear=0;
				sum = 0;
				if(land_mask(y,x) == 0 && invalid_mask(y,x) == 0){
					max_val = samples(y,x,0);
					for(i=-w;i<w+1;i++){
						for(j=-w;j<w+1;j++){
							for(t=0;t<time_size;t++){
								if(land_mask(y+i,x+j) == 0 && invalid_mask(y+i,x+j) == 0){
									if(max_val < samples(y+i,x+j,t)){
										max_val = samples(y+i,x+j,t);
									}

									
									if(std::isnan(clear_samples(y+i,x+j,t))){
										if(!std::isnan(approx_samples(y+i,x+j,t))){										
											sum+= MU_APPROX*approx_samples(y+i,x+j,t);
											count_approx+=1;
										}									
									}
									else{
										sum+=MU_CLEAR*clear_samples(y+i,x+j,t);
										count_clear+=1;
									}
								}
							}
						}
					}
					
					if(count_clear > 0 || count_approx >0 ){
						collated_val = sum/(MU_CLEAR*count_clear +MU_APPROX*count_approx);
						//approx_val = approx_interp(y,x,ind);
						if((max_val - collated_val < T_COLLATED || max_val - collated_val < 0)){// && (approx_val-collated_val < T_COLLATED || approx_val-collated_val < 0 )){
							final_sample(y,x) = collated_val;
							clear(y,x) = count_clear;
							approx(y,x) = count_approx;
						}
					}
				}
			}
		}
}

/*
void
generate_quality_flag(const  Mat1i clear,const Mat1i approx,const Mat1b land_mask,const Mat1b invalid_mask,Mat1i quality_flag)
{
	int qf[SUBSIZE+1][SUBSIZE+1] = {
	 {0, 3, 4, 5, 5, 5, 5, 5},
     {1, 3, 4, 5, 5, 5, 5, 0},
     {2, 3, 4, 5, 5, 5, 0, 0},
     {2, 4, 5, 5, 5, 0, 0, 0},
     {2, 4, 5, 5, 0, 0, 0, 0},
     {2, 4, 5, 0, 0, 0, 0, 0},
     {3, 4, 0, 0, 0, 0, 0, 0},
     {3, 0, 0, 0, 0, 0, 0, 0}
 	};

	int x,y;
	for(y=0;y<HEIGHT;y++){
		for(x=0;x<WIDTH;x++){
			quality_flag(y,x) = qf[approx(y,x)][clear(y,x)];
		}
	}
}
*/
void
subsample(const vector<string> &approx_paths, const vector<string> &original_paths, const vector<string> &smooth_paths, vector<string> &collated_paths,
	      const vector<string> &clear_paths, const Mat1b &land_mask, const Mat1b &invalid_mask)
{
	int i,j;
	int time_size = 7;
	int file_count = 3;
	int mid = time_size/2;
	int dims[3] = {HEIGHT,WIDTH,time_size};
	int approx_dims[3] = {HEIGHT,WIDTH,approx_paths.size()};
	string filename, save_loc, temp_filename;
	string folder_loc = "data/collated";
	int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG+SMOOTH_WINDOW_LAG;

	
	Mat1f bt11_interp_approx(3,approx_dims);
	Mat1f bt12_interp_approx(3,approx_dims);
	Mat1f bt08_interp_approx(3,approx_dims);
	Mat1f bt10_interp_approx(3,approx_dims);
	Mat1f sst_interp_approx(3,approx_dims);
    
    Mat1f bt11_approx(3,dims);
	Mat1f bt12_approx(3,dims);
	Mat1f bt08_approx(3,dims);
	Mat1f bt10_approx(3,dims);
	Mat1f sst_approx(3,dims);

    Mat1f bt11(3,dims);
	Mat1f bt12(3,dims);
	Mat1f bt08(3,dims);
	Mat1f bt10(3,dims);
	Mat1f sst(3,dims);

	Mat1f bt11_clear(3,dims);
	Mat1f bt12_clear(3,dims);
	Mat1f bt08_clear(3,dims);
	Mat1f bt10_clear(3,dims);
	Mat1f sst_clear(3,dims);

	Mat1b clear_masks(3,dims);

	Mat1b bt08_clear_count(HEIGHT,WIDTH);
	Mat1b bt10_clear_count(HEIGHT,WIDTH);
	Mat1b bt11_clear_count(HEIGHT,WIDTH);
	Mat1b bt12_clear_count(HEIGHT,WIDTH);
	Mat1b sst_clear_count(HEIGHT,WIDTH);


    Mat1b bt08_approx_count(HEIGHT,WIDTH);
    Mat1b bt10_approx_count(HEIGHT,WIDTH);
    Mat1b bt11_approx_count(HEIGHT,WIDTH);
    Mat1b bt12_approx_count(HEIGHT,WIDTH);
    Mat1b sst_approx_count(HEIGHT,WIDTH);

	Mat1f bt10_final(HEIGHT,WIDTH);
	Mat1f bt11_final(HEIGHT,WIDTH);
	Mat1f bt12_final(HEIGHT,WIDTH);
	Mat1f bt08_final(HEIGHT,WIDTH);
	Mat1f sst_final(HEIGHT,WIDTH);
    
    /*
	printf("reading approximation files\n");
	for(i=0;i<approx_paths.size();i++){
		readgranule_fullbands(approx_paths[i].c_str(), bt11_interp_approx,bt12_interp_approx,bt08_interp_approx,bt10_interp_approx, sst_interp_approx, i);
	}
	printf("read all approximation files\n");

	printf("starting interpolation of approximation\n");
	interpolate(bt08_interp_approx, land_mask, invalid_mask, time_size, 1000);
	interpolate(bt10_interp_approx, land_mask, invalid_mask, time_size, 1000);
	interpolate(bt11_interp_approx, land_mask, invalid_mask, time_size, 1000);
	interpolate(bt12_interrp_approx, land_mask, invalid_mask, time_size, 1000);
	interpolate(sst_interp_approx, land_mask, invalid_mask, time_size, 1000);
	printf("finished interpolation of approximation\n");
    */

	for(i=mid;i<approx_paths.size() - mid;i++){
    	temp_filename = generate_filename(original_paths[i+original_lag]);
    	if(temp_filename[11] == '0' && temp_filename[12] == '0'){
    		printf("Found time file %s\n",temp_filename.c_str());
			//fill window
			for(j=0;j<time_size;j++){
				

				readgranule_fullbands(original_paths[j+i-mid+original_lag].c_str(), bt11_clear,bt12_clear,bt08_clear,bt10_clear, sst_clear, j);
                readgranule_fullbands(original_paths[j+i-mid+original_lag].c_str(), bt11,bt12,bt08,bt10, sst, j);
                readgranule_fullbands(approx_paths[j+i-mid].c_str(), bt11_approx,bt12_approx,bt08_approx,bt10_approx, sst_approx, j);

				read_mask(clear_paths[j+i-mid+SMOOTH_WINDOW_LAG].c_str(),clear_masks,j);
				apply_mask_slice(clear_masks,bt08_clear,j,true);
				apply_mask_slice(clear_masks,bt10_clear,j,true);
				apply_mask_slice(clear_masks,bt11_clear,j,true);
				apply_mask_slice(clear_masks,bt12_clear,j,true);
				apply_mask_slice(clear_masks,sst_clear,j,true);
	    		printf("loaded file %s\n",original_paths[j+i-mid+original_lag].c_str());
			}

			//generate final
			calc_subsample_3d(bt08_approx, bt08_interp_approx, bt08_clear, bt08, land_mask, invalid_mask, bt08_final, bt08_clear_count, bt08_approx_count, time_size,i, mid);
			printf("processed bt08\n");
			calc_subsample_3d(bt10_approx, bt10_interp_approx,bt10_clear, bt10, land_mask, invalid_mask, bt10_final, bt10_clear_count, bt10_approx_count, time_size,i, mid);
			printf("processed bt10\n");
			calc_subsample_3d(bt11_approx, bt11_interp_approx,bt11_clear, bt11, land_mask, invalid_mask, bt11_final, bt11_clear_count, bt11_approx_count, time_size,i, mid);
			printf("processed bt11\n");
			calc_subsample_3d(bt12_approx, bt12_interp_approx,bt12_clear, bt12, land_mask, invalid_mask, bt12_final, bt12_clear_count, bt12_approx_count, time_size,i, mid);
			printf("processed bt12\n");
			calc_subsample_3d(sst_approx, sst_interp_approx,sst_clear, sst, land_mask, invalid_mask, sst_final, sst_clear_count, sst_approx_count, time_size,i, mid);
			printf("processed sst\n");

			filename = generate_filename(temp_filename);
	    	printf("filename = %s\n",filename.c_str());
	    	file_count+=7;

	    	save_loc = folder_loc+filename;
			save_test_nc_final(bt11_final, bt12_final, bt08_final, bt10_final, sst_final, 
				               bt08_clear_count,bt08_approx_count,bt10_clear_count,bt10_approx_count,bt11_clear_count,bt11_approx_count,
				               bt12_clear_count,bt12_approx_count,sst_clear_count,sst_approx_count,save_loc.c_str());
			collated_paths.push_back(save_loc);
			printf("Generated file: %s\n",save_loc.c_str());
		}
	}
}
