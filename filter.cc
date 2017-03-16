// This function generates assorted cloud masks
// then combines and applies the generated masks
// then performs the second pass 
void 
generate_masks(const Mat1f &bt08, const Mat1f &bt10, const Mat1f &bt11, const Mat1f &bt12, const Mat1f &sst, const Mat1f &reference, const Mat1b &land_mask, 
                const Mat1b &invalid_mask, const Mat1b &ice_mask, const Mat1b &border_mask, Mat1b &bt_mask,int cur_position, string filename)
{

    int x,y;
    Mat1f bt11_clear(HEIGHT,WIDTH);
    Mat1f diffs08(HEIGHT,WIDTH);
    Mat1f diffs11(HEIGHT,WIDTH);
    Mat1f bt_ratio(HEIGHT,WIDTH);

    Mat1f tl0_mask(HEIGHT,WIDTH);
    
    //first copy bt11[:,:,ind] into bt11_clear[:,:,ind]
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            bt11_clear(y,x) = sst(y,x,cur_position);
            tl0_mask(y,x) = sst(y,x,cur_position);
            //cold_mask(y,x) = bt11(y,x,cur_position);
            //ratio_mask(y,x) = bt11(y,x,cur_position);
            //nn_mask(y,x) = bt11(y,x,cur_position);
            //diag_mask(y,x) = bt11(y,x,cur_position);  
            //t_mask(y,x) = bt11(y,x,cur_position);           
        }
    }



    apply_l2p_flags(ice_mask, land_mask, invalid_mask, bt11_clear, cur_position,true);
    //printf("starting histogram\n");
    //filter_histogram(bt08, bt10,bt11,bt12, lats, bt11_clear, cur_position);
    //printf("finished histogram\n");

    
    compute_cold_mask(bt08,bt10,bt11_clear,cur_position);
    printf("finished coldmask\n");
    //SAVENC(cold_mask);      
    //compute_cold_mask(bt08,bt10,cold_mask,cur_position);

    compute_dtmask(reference,sst,bt11_clear,cur_position, T_DT);
    //compute_threshmask(dt_analysis,dt_mask,cur_position, T_DT, true);

    calculate_bt_ratio(bt08, bt11,bt12, bt_ratio, RATIO, cur_position);
    //SAVENC(bt_ratio);
    compute_threshmask_2d(bt_ratio, bt11_clear,cur_position,T_RATIO,true);
    //compute_threshmask_2d(bt_ratio, ratio_mask,cur_position,T_RATIO,true);
    //printf("finished bt11 diffs mask\n");
    //SAVENC(bt11_mask);
    
    
    compute_nnmask(bt11, bt11_clear,cur_position);
    //compute_nnmask(bt11, nn_mask,cur_position);
    printf("finished nnmask\n");
    //SAVENC(nn_mask);  
    
    //compute_avgmask(bt11,bt11_clear);

    compute_diagmask(bt11, bt11_clear,cur_position);
    //compute_diagmask(bt11, diag_mask,cur_position);
    printf("finished diagmask\n");
    //SAVENC(diag_mask);
    compute_tl0(sst, border_mask, bt11_clear,cur_position);
    //compute_tl0(sst, border_mask, tl0_mask,cur_position);
    
    //compute_eigenvals(bt08,bt10,bt11,bt12, border_mask,bt11_clear,cur_position);
    printf("finished eigen mask\n");
    
    

    /*
    generate_mask(cold_mask, bt_mask);
    string save_loc = "data/cold_mask" + filename;
    save_test_nc_float(bt_mask,save_loc.c_str());

    generate_mask(dt_mask, bt_mask);
    save_loc = "data/dt_mask" + filename;
    save_test_nc_float(bt_mask,save_loc.c_str());

    generate_mask(ratio_mask, bt_mask);
    save_loc = "data/ratio_mask" + filename;
    save_test_nc_float(bt_mask,save_loc.c_str());

    generate_mask(nn_mask, bt_mask);
    save_loc = "data/nn_mask" + filename;
    save_test_nc_float(bt_mask,save_loc.c_str());

    generate_mask(diag_mask, bt_mask);
    save_loc = "data/diag_mask" + filename;
    save_test_nc_float(bt_mask,save_loc.c_str());
    */

    //generate_mask(tl0_mask, bt_mask);
    //string save_loc = "data/tl0_mask" + filename;
    //save_test_nc_float(bt_mask,save_loc.c_str());

    generate_mask(bt11_clear, bt_mask);
    
    //SAVENC(bt11_clear);
    /*
    printf("starting 2nd pass\n");
    clear_samples_compare(bt11,bt11_clear, land_mask, invalid_mask,bt11_masked,dims, inds, PASS2);
    printf("finished 2nd pass\n");
    */
}




void
filter_clouds(const vector<string> &paths, const Mat1b &land_mask, const Mat1b &invalid_mask, const Mat1b &border_mask, 
    vector<string> &clear_paths, const string ref_file){

    int i=0;
    int j;
    int sample_size = paths.size();
    //int sample_size = 19;
    //int sample_size = 27;
    vector<string> curPaths;      
    int time_ind =0 ;
    int new_granule_ind = 0;
    int current_mid = FILTER_TIME_SIZE/2;
    int new_mask_ind = current_mid+1;
    int file_count = current_mid;
    int x,y,t;
    int dims[3] = {HEIGHT,WIDTH,FILTER_TIME_SIZE};
    int inds[FILTER_TIME_SIZE];
    // Cur controls where the mask is read from

    string filename;
    string folder_loc = "data/clear";
    string save_loc;

    Mat1f bt11(3,dims);
    Mat1f bt12(3,dims);
    Mat1f bt08(3,dims);
    Mat1f bt10(3,dims);
    Mat1f sst(3,dims);
    Mat1f reference(HEIGHT,WIDTH);
    Mat1b bt_mask(HEIGHT, WIDTH);
    Mat1f temp(HEIGHT, WIDTH);
    Mat1f tempsst(HEIGHT, WIDTH);
    Mat1f bt11_clear();
    Mat1b ice_masks(3,dims);

    read_reference(ref_file,reference,"sst_reynolds");

    for(j=0;j<FILTER_TIME_SIZE;j++){

        readgranule(paths[time_ind], bt11,bt12,bt08,bt10, sst, j);
        get_icemask(paths[time_ind],ice_masks, j);
        apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt11, j,false);
        apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt12, j,false);
        apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt10, j,false);
        apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt08, j,false);
        apply_l2p_flags(ice_masks, land_mask, invalid_mask, sst, j,false);
        time_ind++;
    }


    //SAVENC(dt_analysis);
    filename = generate_filename(paths[file_count]);
    generate_masks(bt08,bt10,bt11,bt12, sst,reference,land_mask,invalid_mask, ice_masks, border_mask, bt_mask,current_mid,filename);
    
    
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
            readgranule(paths[time_ind], bt11,bt12,bt08,bt10, sst, new_granule_ind);
            printf("read new granule\n");
            get_icemask(paths[time_ind], ice_masks,new_granule_ind);
            printf("read new ice mask\n");
            apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt11, new_granule_ind,false);
            apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt12, new_granule_ind,false);
            apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt10, new_granule_ind,false);
            apply_l2p_flags(ice_masks, land_mask, invalid_mask, bt08, new_granule_ind,false);
            apply_l2p_flags(ice_masks, land_mask, invalid_mask, sst, new_granule_ind,false);
            new_granule_ind = (new_granule_ind+1) % FILTER_TIME_SIZE;
            printf("starting mask generation\n");
            filename = generate_filename(paths[file_count]);
            generate_masks(bt08,bt10,bt11,bt12,sst, reference,land_mask,invalid_mask, ice_masks,border_mask, bt_mask,new_mask_ind,filename);
            new_mask_ind = (new_mask_ind+1) % FILTER_TIME_SIZE;  
            //clear_samples_compare(bt11, bt11_clear, land_mask, invalid_mask, ice_masks, bt11_masked,current_mid);
            current_mid =(current_mid + 1) % FILTER_TIME_SIZE; 
            
            file_count++;         
        }
        time_ind++;i++;
    }
    
}
       
    