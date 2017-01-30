typedef struct stat
{
    float sum=0;
    int count=0;
}stat;

void
smoothing3D_std(const Mat1f &clear, const Mat1b &land_mask, const Mat1b &invalid_mask, Mat1f &sst_3dsmooth, int *inds)
{
	int x,y,z,i,j,k;
	int t_lag = 9;
	int w = 7;

	int t_len = 2*t_lag+1;
	int w_len = 2*w+1;
	int threshold;
	double left_sum, sum, final_sum, val, avg, stdev, check, count, left_count, final_count;
	double sqsum, left_sqsum, sqavg;

	deque<double> nonNan_values;

	threshold = PR_CLEAR * (2*t_lag+1)*(2*w+1)*(2*w+1);
	//threshold =1;
	//printf("%d \n",threshold);
	z = t_lag;
	printf("starting loop\n");

	//calculate stats in time dimension
	for(y=w;y<HEIGHT-w+1;y++){
		x=0;
		sum = 0;
		count  = 0;
		left_sum = 0;
		left_count = 0;
		sqsum = 0;
		left_sqsum =0;
		nonNan_values.clear();
		//find first first valid pixel
		while(land_mask(y,x) !=0 && invalid_mask(y,x) != 0){
			x++;
		}
		
        // calc rest of window sum and count		
		for(i=-w;i<w+1;i++){
			for(j=-w;j<w+1;j++){
				for(k=-t_lag;k<t_lag+1;k++){
					val = clear(y+i,x+j,inds[z+k]);
					if(!std::isnan(val)){
						count+=1;
						sum += val;
						sqsum += val * val;
						nonNan_values.push_back(val);
					}
				}
			}
		}

		//calculate average and determine what to use for smooth
		final_count =0;
		final_sum = 0;
		if(count >= threshold){
			avg = sum/count;
			sqavg = sqsum/count;
			stdev = sqrt(sqavg - avg*avg);
			for(i=0;i<nonNan_values.size();i++){
				check = fabs((nonNan_values[i] - avg)/stdev);
				if(check < 2){
					final_sum += nonNan_values[i];
					final_count++;
				}
			}
		}

		if(final_count >0){

			sst_3dsmooth(y,x) = final_sum /final_count;
		}

		// now repeat for all other valid x values
		for(x=x+1;x<WIDTH;x++){
			//remove all values of previous left
			sum -= left_sum;
			count -= left_count;
			sqsum -= left_sqsum;

			for(i=0;i<left_count;i++){
				nonNan_values.pop_front();
			}

			//calculate new left
			left_sum = 0;
			left_count = 0;
			left_sqsum = 0;
			j = -w;
			for(i=-w;i<w+1;i++){
				for(k=-t_lag;k<t_lag+1;k++){
					val = clear(y+i,x+j,inds[z+k]);
					if(!std::isnan(val)){
						left_count+=1;
						left_sum += val;
						left_sqsum += val*val;
					}
				}
			}

			//calculate new right
			j = w;
			for(i=-w;i<w+1;i++){
				for(k=-t_lag;k<t_lag+1;k++){
					val = clear(y+i,x+j,inds[z+k]);
					if(!std::isnan(val)){
						count+=1;
						sum += val;
						sqsum += val*val;
						nonNan_values.push_back(val);
					}
				}
			}
			//printf("count = %f sum = %f\n", count, sum);
			if(land_mask(y,x)==0 && invalid_mask(y,x) ==0){
				// calculate this pixels value at 3dsmooth(y,x)
				final_count =0;
				final_sum = 0;
				if(count >= threshold){
					avg = sum/count;
					sqavg = sqsum/count;
					stdev = sqrt(sqavg - (avg*avg));

					for(i=0;i<nonNan_values.size();i++){
						check = fabs((nonNan_values[i] - avg)/stdev);
						if(check < 2){
							final_sum += nonNan_values[i];
							final_count++;
						}
					}
				}
				if(final_count >0){
					sst_3dsmooth(y,x) = final_sum /final_count;
				}
			}
		}

	}

}


void
smoothing3D_std(const Mat1f &clear, const Mat1b &land_mask, const Mat1b &invalid_mask, Mat1f &sst_3dsmooth, int *inds)
{
	int x,y,z,i,j,k;
	int t_lag = 9;
	int w = 7;

	int t_len = 2*t_lag+1;
	int w_len = 2*w+1;
	int threshold;
	double left_sum, sum,  val, avg,  count, left_count;


	threshold = PR_CLEAR * (2*t_lag+1)*(2*w+1)*(2*w+1);
	//threshold =1;
	//printf("%d \n",threshold);
	z = t_lag;
	printf("starting loop\n");

	//calculate stats in time dimension
	for(y=w;y<HEIGHT-w+1;y++){
		x=0;
		sum = 0;
		count  = 0;
		left_sum = 0;
		left_count = 0;
		//find first first valid pixel
		while(land_mask(y,x) !=0 && invalid_mask(y,x) != 0){
			x++;
		}
		
        // calc rest of window sum and count		
		for(i=-w;i<w+1;i++){
			for(j=-w;j<w+1;j++){
				for(k=-t_lag;k<t_lag+1;k++){
					val = clear(y+i,x+j,inds[z+k]);
					if(!std::isnan(val)){
						count+=1;
						sum += val;
					}
				}
			}
		}

		if(count > threshold){

			sst_3dsmooth(y,x) = sum /count;
		}

		// now repeat for all other valid x values
		for(x=x+1;x<WIDTH;x++){
			//remove all values of previous left
			sum -= left_sum;
			count -= left_count;

			//calculate new left
			left_sum = 0;
			left_count = 0;
			j = -w;
			for(i=-w;i<w+1;i++){
				for(k=-t_lag;k<t_lag+1;k++){
					val = clear(y+i,x+j,inds[z+k]);
					if(!std::isnan(val)){
						left_count+=1;
						left_sum += val;
					}
				}
			}

			//calculate new right
			j = w;
			for(i=-w;i<w+1;i++){
				for(k=-t_lag;k<t_lag+1;k++){
					val = clear(y+i,x+j,inds[z+k]);
					if(!std::isnan(val)){
						count+=1;
						sum += val;
					}
				}
			}
			//printf("count = %f sum = %f\n", count, sum);
			if(land_mask(y,x)==0 && invalid_mask(y,x) ==0){
				// calculate this pixels value at 3dsmooth(y,x)
				if(count>threshold){
					sst_3dsmooth(y,x) = sum /count;
				}
			}
		}
	}
}



void
smooth_samples(const vector<string> paths, const Mat1b &land_mask, const Mat1b &invalid_mask,vector<string> &smooth_paths){
	
	int j;
	int i=0;
    int time_size = 19;
    int time_ind = time_size;
    int start_ind = 0;
    int dims[3] = {HEIGHT,WIDTH,time_size};
    int inds[time_size];
    //int sample_size = paths.size();
    int sample_size = 19;

    Mat1f clear_samples(3,dims);
  	Mat1f bt11_3dsmooth(HEIGHT, WIDTH);

  	compute_indicies(inds, start_ind, time_size);

    for(j=0;j<time_size;j++){
    	readgranules_oneband(paths[j],clear_samples, j);
    	
    }
    printf("read granules\n");
    while(time_ind < sample_size+1){     
 
        printf("starting 3d smoothing \n");
       
        smoothing3D(clear_samples, land_mask, invalid_mask, bt11_3dsmooth, inds);
        printf("finished 3d smoothing\n");
               
        
        string filename = "data/smooth/smooth_samples_small" + convert_int_to_string(i) +".nc";
        save_test_nc_float(bt11_3dsmooth,filename.c_str());        
        smooth_paths.push_back(filename.c_str());
        printf("Generated File %s\n",filename.c_str());
        
        if(time_ind<sample_size){
	        readgranules_oneband(paths[time_ind],clear_samples,start_ind);
	        printf("read granule %d\n",time_ind);        
	        start_ind = (start_ind +1) % time_size;        
	        printf("starting compute_indicies\n");
	        compute_indicies(inds,start_ind,time_size);
	        printf("finished compute_indicies\n");
	        
	      }
	    time_ind++; i++;
    }
        
}
 

