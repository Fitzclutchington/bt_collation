// This function generates assorted cloud masks
// then combines and applies the generated masks
// then performs the second pass 
void 
generate_masks(const Mat1f &bt08, const Mat1f &bt10, const Mat1f &bt11, const Mat1f &bt12, const Mat1f &dt_analysis, const Mat1b &land_mask, 
                const Mat1b &invalid_mask, const Mat1b &ice_mask, const Mat1b &border_mask, Mat1b &bt_mask,int cur_position)
{

    int x,y;
    Mat1f bt11_clear(HEIGHT,WIDTH);
    Mat1f diffs08(HEIGHT,WIDTH);
    Mat1f diffs11(HEIGHT,WIDTH);
    Mat1f bt_ratio(HEIGHT,WIDTH);

    //first copy bt11[:,:,ind] into bt11_clear[:,:,ind]
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            bt11_clear(y,x) = bt11(y,x,cur_position);            
        }
    }



    apply_l2p_flags(ice_mask, land_mask, invalid_mask, bt11_clear, cur_position);
    //printf("starting histogram\n");
    //filter_histogram(bt08, bt10,bt11,bt12, lats, bt11_clear, cur_position);
    //printf("finished histogram\n");

    
    compute_cold_mask(bt08,bt10,bt11_clear,cur_position);
    printf("finished coldmask\n");
    //SAVENC(cold_mask);      
    
    compute_threshmask(dt_analysis,bt11_clear,cur_position, T_DT, true);
    calculate_bt_ratio(bt08, bt11,bt12, bt_ratio, RATIO, cur_position);
    //SAVENC(bt_ratio);
    compute_threshmask_2d(bt_ratio, bt11_clear,cur_position,T_RATIO,true);
    //printf("finished bt11 diffs mask\n");
    //SAVENC(bt11_mask);
    
    
    compute_nnmask(bt11, bt11_clear,cur_position);
    printf("finished nnmask\n");
    //SAVENC(nn_mask);  
    
    
    compute_diagmask(bt11, bt11_clear,cur_position);
    printf("finished diagmask\n");
    //SAVENC(diag_mask);
    
    compute_eigenvals(bt08,bt10,bt11,bt12, border_mask,bt11_clear,cur_position);
    printf("finished eigen mask\n");
    
    generate_mask(bt11_clear, bt_mask);
    //SAVENC(bt11_clear);
    /*
    printf("starting 2nd pass\n");
    clear_samples_compare(bt11,bt11_clear, land_mask, invalid_mask,bt11_masked,dims, inds, PASS2);
    printf("finished 2nd pass\n");
    */
}




void
filter_clouds(const vector<string> &paths, const Mat1b &land_mask, const Mat1b &invalid_mask, const Mat1b &border_mask, vector<string> &clear_paths){

    int i=0;
    int j;
    int sample_size = paths.size();
    //int sample_size = 19;
    //int sample_size = 27;
    vector<string> curPaths;      
    int time_ind =0 ;
    int new_granule_ind = 0;
    int new_mask_ind = FILTER_TIME_SIZE -1;
    int current_mid = new_mask_ind/2;
    int file_count = 1;
    int dims[3] = {HEIGHT,WIDTH,FILTER_TIME_SIZE};
    // Cur controls where the mask is read from

    string filename;
    string folder_loc = "data/clear";
    string save_loc;

    Mat1f bt11(3,dims);
    Mat1f bt12(3,dims);
    Mat1f bt08(3,dims);
    Mat1f bt10(3,dims);
    Mat1f dt_analysis(3,dims);
    Mat1b bt_mask(HEIGHT, WIDTH);
    Mat1f temp(HEIGHT, WIDTH);
    Mat1f tempsst(HEIGHT, WIDTH);
    Mat1f bt11_clear();
    Mat1b ice_masks(3,dims);

    for(j=0;j<FILTER_TIME_SIZE;j++){

        readgranule(paths[time_ind], bt11,bt12,bt08,bt10, dt_analysis, j);
        get_icemask(paths[time_ind],ice_masks, j);
        time_ind++;
    }

    //SAVENC(dt_analysis);
    for(j=FILTER_WINDOW_LAG;j<FILTER_TIME_SIZE-FILTER_WINDOW_LAG;j++){
        generate_masks(bt08,bt10,bt11,bt12, dt_analysis,land_mask,invalid_mask, ice_masks, border_mask, bt_mask,j);
    }
    
    //bt11_masked.setTo(NAN);
    //clear_samples_compare(bt11, bt11_clear, land_mask,invalid_mask, ice_masks, bt11_masked,current_mid);
    //

    filename = generate_filename(paths[file_count]);
    printf("filename = %s\n",filename);
    file_count++;

    while(time_ind < sample_size+1){
        // can be sped up by maintaining a dequeue of all data
        // pop and push at next iteration
        
        //printf("starting cloud mask \n");
        //clear_samples_compare(bt11,bt11_clear, land_mask, invalid_mask,bt11_masked,inds, PASS2);           
        //printf("finished cloud mask\n");
        

        save_loc = folder_loc+ filename;       
        

        // make this use byte type
        //SAVENC(bt_mask);
        save_test_nc_float(bt_mask,save_loc.c_str());

        clear_paths.push_back(save_loc.c_str());
        printf("Generated File %s\n",save_loc.c_str());
        
        if(time_ind < sample_size){
            printf("continuing loop\n");
            readgranule(paths[time_ind], bt11,bt12,bt08,bt10, dt_analysis, new_granule_ind);
            printf("read new granule\n");
            get_icemask(paths[time_ind], ice_masks,new_granule_ind);
            printf("read new ice mask\n");
            new_granule_ind = (new_granule_ind+1) % FILTER_TIME_SIZE;
            printf("starting mask generation\n");
            generate_masks(bt08,bt10,bt11,bt12, dt_analysis,land_mask,invalid_mask, ice_masks,border_mask, bt_mask,new_mask_ind);
            new_mask_ind = (new_mask_ind+1) % FILTER_TIME_SIZE;  
            //clear_samples_compare(bt11, bt11_clear, land_mask, invalid_mask, ice_masks, bt11_masked,current_mid);
            current_mid =(current_mid + 1) % FILTER_TIME_SIZE; 
            filename = generate_filename(paths[file_count]);
            file_count++;         
        }
        time_ind++;i++;
    }
    
}
       
    
