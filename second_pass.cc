// bt11_clear is an MxNxT matrix containing the clear samples
// bt11_masked is an MxN matrix which is the slice being returned containing the final mask
void
clear_samples_compare(const Mat1f &bt11, const Mat1f &bt11_clear, const Mat1b &land_mask, const Mat1b &invalid_mask,
                      Mat1b &bt11_masked,int cur_mid, int iter)
{
    int x,y,t,time;
    int t_lag0 = 6;
    int t_lag1 = 2;
    //int t_lag0=3;           //(one hour padding on each side)
    //int t_lag1=1;           
    int w1=2;
    int w0 = 0;
    int t_lag2=2*t_lag0;    //(2 hours padding on each side)
    int w2=2;               // spatial window 2*w1+1 
    
    string filename;
    string num;

    Mat1f TL1(HEIGHT,WIDTH);
    Mat1f TL2(HEIGHT,WIDTH);
    Mat1f TL0(HEIGHT,WIDTH);
    Mat1f t_lag0_count(HEIGHT,WIDTH);

    TL1.setTo(NAN);
    TL2.setTo(NAN);
    TL0.setTo(NAN);

    t_lag0_count.setTo(0);

    int window0[3] = {w0,w0,t_lag1};
    int window1[3] = {w1,w1,t_lag1};
    int window2[3] = {w2,w2,t_lag2};
    
                    //calculate NANmean in data window
    
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            for(t=-t_lag0;t<t_lag0+1;t++){
                time =(((cur_mid + t) % SECOND_PASS_SIZE + SECOND_PASS_SIZE)%SECOND_PASS_SIZE);
                //printf("time = %d\n",time);
                if(!std::isnan(bt11_clear(y,x,time))){
                    t_lag0_count(y,x)+=1;
                }
            }
        }
    }

    windowed_nanmean_2nd_pass(bt11, t_lag0_count, land_mask,invalid_mask, TL1, window1, cur_mid);

    printf("finished first window\n");

    windowed_nanmean_2nd_pass(bt11_clear, t_lag0_count, land_mask,invalid_mask, TL2, window2, cur_mid);

    printf("finished second window\n");

    windowed_nanmean_2nd_pass(bt11, t_lag0_count, land_mask,invalid_mask, TL0, window0, cur_mid);

    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(!std::isnan(TL1(y,x)) && !std::isnan(TL2(y,x))){
                if( bt11(y,x,cur_mid) > TL0(y,x) - PASS2 && ((TL2(y,x)-TL1(y,x) < PASS2 && !std::isnan(bt11_clear(y,x,cur_mid))))){ //|| bt11(y,x,cur_mid) - TL2(y,x) > PASS2 )){
                    bt11_masked(y,x) = 255;
                }
                //if( !std::isnan(bt11_clear(y,x,cur_mid))){
                //    bt11_masked(y,x) = 255;
                //}
            }
        }
    }
    
    /*
    filename = "data/TL1/TL1_" + convert_int_to_string(iter) +".nc";
    save_test_nc_float(TL1,filename.c_str());
    filename = "data/TL2/TL2_" + convert_int_to_string(iter) + ".nc";
    save_test_nc_float(TL2,filename.c_str());
    filename = "data/TL0/TL0_" + convert_int_to_string(iter) + ".nc";
    save_test_nc_float(TL0,filename.c_str());
    */
}    

/* function deals with opening files, passes matrix to second passfunction, and 
   returns new mask to file */
void
second_pass(const vector<string> &clear_paths, const vector<string> &original_paths, vector<string> &second_pass_paths, Mat1b &land_mask, Mat1b &invalid_mask)
{
	int j;
	int time_ind = 0;
	int sample_size = clear_paths.size();
	int time_size = SECOND_PASS_SIZE;
	int dims[3] = {HEIGHT, WIDTH, time_size};
	int new_granule_ind = 0;
    int cur_mid = time_size /2;
    int file_count = cur_mid;
    int iter = 0;

	string folder_loc = "data/pass2";
	string filename, save_loc;
	string variable_name = "sea_surface_temperature";
	string original_filenames[25];
	string clear_filenames[25];

	Mat1b clear_mask(HEIGHT,WIDTH);
	Mat1b clear2(HEIGHT,WIDTH);
	Mat1f sst_samples(3,dims);
    Mat1f clear_samples(3,dims);
   

    clear2.setTo(0);

	for(j=0;j<time_size;j++){
		// read granule one band
    	read_mask(clear_paths[j],clear_mask,-1);
        //read_acspo(original_paths[j+FILTER_WINDOW_LAG],clear_mask,-1);
      
    	readgranule_oneband(original_paths[j+FILTER_WINDOW_LAG],sst_samples,j,variable_name);
    	readgranule_oneband(original_paths[j+FILTER_WINDOW_LAG],clear_samples,j,variable_name);
    	apply_mask_slice(clear_mask, clear_samples, j,false);
    	//get_icemask(original_paths[time_ind+FILTER_WINDOW_LAG],ice_masks, j);
    	original_filenames[j] = original_paths[j+FILTER_WINDOW_LAG];
    	clear_filenames[j] = clear_paths[j];
    	time_ind++;
    	//printf("read clear at %d\n",time_ind);
    	//printf("original file = %s\n",original_paths[j+FILTER_WINDOW_LAG].c_str());
    	//printf("clear file = %s\n",clear_paths[j].c_str());
    }
    
    clear_samples_compare(sst_samples, clear_samples, land_mask, invalid_mask, clear2,cur_mid, iter);
    iter++;
 	filename = generate_filename(clear_paths[file_count]);
    file_count++;

    while(time_ind < sample_size+1){
        // can be sped up by maintaining a dequeue of all data
        // pop and push at next iteration
        
        //printf("starting cloud mask \n");
        //clear_samples_compare(bt11,bt11_clear, land_mask, invalid_mask,bt11_masked,inds, PASS2);           
        //printf("finished cloud mask\n");
        
    	//printf("current mid = %d\n", cur_mid);
        save_loc = folder_loc+ filename;       
        

        // make this use byte type
        //SAVENC(bt_mask);
        save_test_nc_float(clear2,save_loc.c_str());


        second_pass_paths.push_back(save_loc.c_str());
        printf("Generated File %s\n",save_loc.c_str());
        printf("finished iteration %d\n",iter);
        if(time_ind < sample_size){
            read_mask(clear_paths[time_ind],clear_mask,-1);
            //read_acspo(original_paths[time_ind+FILTER_WINDOW_LAG],clear_mask,-1);
            //printf("original path name = %s\n", original_paths[time_ind + FILTER_WINDOW_LAG].c_str());
            //printf("clear path name = %s\n", clear_paths[time_ind].c_str());
    		readgranule_oneband(original_paths[time_ind+FILTER_WINDOW_LAG],sst_samples,new_granule_ind,variable_name);
    		readgranule_oneband(original_paths[time_ind+FILTER_WINDOW_LAG],clear_samples,new_granule_ind,variable_name);
    		apply_mask_slice(clear_mask, clear_samples, new_granule_ind,false);
    		//get_icemask(original_paths[time_ind+FILTER_WINDOW_LAG],ice_masks, new_granule_ind);

    		original_filenames[new_granule_ind] = original_paths[time_ind+FILTER_WINDOW_LAG];
    		clear_filenames[new_granule_ind] = clear_paths[time_ind];

            new_granule_ind = (new_granule_ind+1) % SECOND_PASS_SIZE;
            cur_mid = (cur_mid + 1) % SECOND_PASS_SIZE;
            printf("starting mask generation\n");
            clear2.setTo(0);
            clear_samples_compare(sst_samples, clear_samples, land_mask, invalid_mask, clear2,cur_mid, iter);             
            filename = generate_filename(clear_paths[file_count]);
            file_count++; iter++;       
        }
        time_ind++;
    }
}
