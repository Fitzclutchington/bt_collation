
void
smoothing3D(const Mat1f &bt08, const Mat1f &bt10, const Mat1f &bt11, const Mat1f &bt12, const Mat1f &sst, const Mat1b &land_mask, const Mat1b &invalid_mask, 
	const Mat1b &ice_masks, Mat1f &bt08_3dsmooth, Mat1f &bt10_3dsmooth, Mat1f &bt11_3dsmooth, Mat1f &bt12_3dsmooth, Mat1f &sst_3dsmooth, int *inds)
{
	int x,y,t,i,j,k;
	int t_lag = 9;
	int w = 10;

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

	windowed_nanmean_3d(bt08, land_mask, invalid_mask, ice_masks, bt08_3dsmooth, window, inds, threshold);
	windowed_nanmean_3d(bt10, land_mask, invalid_mask, ice_masks, bt10_3dsmooth, window, inds, threshold);
	windowed_nanmean_3d(bt11, land_mask, invalid_mask, ice_masks, bt11_3dsmooth, window, inds, threshold);
	windowed_nanmean_3d(bt12, land_mask, invalid_mask, ice_masks, bt12_3dsmooth, window, inds, threshold);
	windowed_nanmean_3d(sst, land_mask, invalid_mask, ice_masks, sst_3dsmooth, window, inds, threshold);
}



void
smooth_samples(const vector<string> clear_paths, const vector<string> original_paths,
			   const Mat1b &land_mask, const Mat1b &invalid_mask,vector<string> &smooth_paths){
	
	int j;
	int i=0;
    int time_size = SMOOTH_TIME_WINDOW;
    int time_ind = time_size;
    int file_count = (time_size -1) /2;
    int start_ind = 0;
    int dims[3] = {HEIGHT,WIDTH,time_size};
    int inds[time_size];
    int sample_size = clear_paths.size();
    string folder_loc = "data/smooth";
 	string filename, save_loc;

 	Mat1f bt08(3,dims);
    Mat1f bt10(3,dims);
 	Mat1f bt11(3,dims);
    Mat1f bt12(3,dims);
    Mat1f sst(3,dims); 

    Mat1b clear_mask(HEIGHT,WIDTH);

  	Mat1f bt08_3dsmooth(HEIGHT, WIDTH);
  	Mat1f bt10_3dsmooth(HEIGHT, WIDTH);
  	Mat1f bt11_3dsmooth(HEIGHT, WIDTH);
  	Mat1f bt12_3dsmooth(HEIGHT, WIDTH);
  	Mat1f sst_3dsmooth(HEIGHT, WIDTH);

  	Mat1b ice_masks(3,dims);

  	compute_indicies(inds, start_ind, time_size);

    for(j=0;j<time_size;j++){
    	
    	readgranule_fullbands(original_paths[j+FILTER_WINDOW_LAG + SECOND_PASS_LAG], bt11,bt12,bt08,bt10, sst, j);
    	printf("read brightness temps\n");
    	read_mask(clear_paths[j],clear_mask,-1);
    	printf("read mask\n");
    	apply_mask(clear_mask,bt08,bt10,bt11,bt12, sst, j);
    	printf("applied mask\n");
    	get_icemask(original_paths[j], ice_masks, j);
    	printf("got ice mask mask\n");    	
    }
    //SAVENC(bt08);
    printf("read granules\n");
    while(time_ind < sample_size+1){     
 
        printf("starting 3d smoothing \n");
        bt08_3dsmooth.setTo(NAN);
        bt10_3dsmooth.setTo(NAN);
        bt11_3dsmooth.setTo(NAN);
        bt12_3dsmooth.setTo(NAN);
        sst_3dsmooth.setTo(NAN);
        smoothing3D(bt08, bt10, bt11, bt12, sst,land_mask, invalid_mask, ice_masks, bt08_3dsmooth, bt10_3dsmooth, bt11_3dsmooth, bt12_3dsmooth, sst_3dsmooth, inds);
        printf("finished 3d smoothing\n");
               
        filename = generate_filename(clear_paths[file_count]);
        save_loc = folder_loc+ filename; 
    	file_count++;

        save_test_nc_fullbands(bt08_3dsmooth, bt10_3dsmooth, bt11_3dsmooth, bt12_3dsmooth, sst_3dsmooth, save_loc.c_str());        
        smooth_paths.push_back(save_loc.c_str());
        printf("Generated File %s\n",save_loc.c_str());
        
        if(time_ind<sample_size){
	        readgranule_fullbands(original_paths[time_ind+FILTER_WINDOW_LAG+SECOND_PASS_LAG], bt11,bt12,bt08,bt10, sst, start_ind);
	        read_mask(clear_paths[time_ind],clear_mask,-1);
    		apply_mask(clear_mask,bt08,bt10,bt11,bt12, sst, time_ind);
	        get_icemask(original_paths[time_ind], ice_masks, start_ind);
	        printf("read granule %d\n",time_ind);        
	        start_ind = (start_ind +1) % time_size;        
	        printf("starting compute_indicies\n");
	        compute_indicies(inds,start_ind,time_size);
	        printf("finished compute_indicies\n");
	        
	      }
	    time_ind++; i++;
    }
        
}