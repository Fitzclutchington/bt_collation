
void
smoothing3D(const Mat1f &clear, const Mat1b &land_mask, const Mat1b &invalid_mask, Mat1f &sst_3dsmooth, int *inds)
{
	int x,y,t,i,j,k;
	int t_lag = 9;
	int w = 7;

	int window[3] = {w,w,t_lag};

	int t_len = 2*t_lag+1;
	int w_len = 2*w+1;
	int threshold;
	double left_sum, sum,  val, avg,  count, left_count;

	Mat1f time_sum(HEIGHT,WIDTH);
	Mat1f time_count(HEIGHT,WIDTH);

	threshold = PR_CLEAR * (2*t_lag+1)*(2*w+1)*(2*w+1);
	//threshold =1;
	//printf("%d \n",threshold);
	t = t_lag;
	printf("starting loop\n");

	windowed_nanmean_3d(clear, land_mask, invalid_mask, sst_3dsmooth, window, inds, threshold);
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
    int sample_size = paths.size();

 

    Mat1f clear_samples(3,dims);
  	Mat1f bt11_3dsmooth(HEIGHT, WIDTH);

  	compute_indicies(inds, start_ind, time_size);

    for(j=0;j<time_size;j++){
    	readgranules_oneband(paths[j],clear_samples, j);
    	
    }
    printf("read granules\n");
    while(time_ind < sample_size+1){     
 
        printf("starting 3d smoothing \n");
        bt11_3dsmooth.setTo(NAN);
        smoothing3D(clear_samples, land_mask, invalid_mask, bt11_3dsmooth, inds);
        printf("finished 3d smoothing\n");
               
        
        string filename = "data/smooth/smooth_samples" + convert_int_to_string(i) +".nc";
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
 

