// This function generates assorted cloud masks
// then combines and applies the generated masks
// then performs the second pass 
void 
generate_masks(const Mat1f &bt08, const Mat1f &bt10, const Mat1f &bt11, const Mat1f &bt12, const Mat1f &sst, 
               const Mat1f &reference, const Mat1b &border_mask, Mat1b &bt_mask,int cur_position)
{

    int x,y;
    Mat1f bt11_clear(HEIGHT,WIDTH);
    Mat1f bt_ratio(HEIGHT,WIDTH);
    
    //first copy bt11[:,:,ind] into bt11_clear[:,:,ind]
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            bt11_clear(y,x) = sst(y,x,cur_position);          
        }
    }

    compute_dtmask(reference,sst,bt11_clear,cur_position, T_DT);
    printf("finished dt dt_analysis\n");

    calculate_bt_ratio(bt08, bt11,bt12, bt_ratio, cur_position);
    compute_threshmask_2d(bt_ratio, bt11_clear,T_RATIO,true);    
    printf("finished 2 bt mask\n");

    compute_nnmask(bt11, bt11_clear,cur_position);
    printf("finished nnmask\n");  
    
    compute_cold_mask(bt08,bt10,bt11_clear,cur_position);
    printf("finished 2 bt Mask\n"); 

    compute_diagmask(bt11, bt11_clear,cur_position);
    printf("finished diagmask\n");
    
    compute_tl0(sst, border_mask, bt11_clear);
    printf("finished eigen mask\n");  

    generate_mask(bt11_clear, bt_mask); 

}




void
filter_clouds(const vector<string> &paths, const Mat1b &l2p_mask, const Mat1b &border_mask, 
              vector<string> &clear_paths, const string ref_file){

    int i=0;
    int j;
    int sample_size = paths.size();
    vector<string> curPaths;      
    int time_ind =0 ;
    int new_granule_ind = 0;
    int current_mid = FILTER_TIME_SIZE/2;
    int new_mask_ind = current_mid+1;
    int file_count = current_mid;
    int dims[3] = {HEIGHT,WIDTH,FILTER_TIME_SIZE};
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

    read_reference(ref_file,reference,"sst_reynolds");

    // Fill data arrays with inital brightness_temperature/sst data
    for(j=0;j<FILTER_TIME_SIZE;j++){
        readgranule(paths[time_ind], bt11,bt12,bt08,bt10, sst, j);
        apply_l2p_flags( l2p_mask, bt11, j,false);
        apply_l2p_flags( l2p_mask, bt12, j,false);
        apply_l2p_flags( l2p_mask, bt10, j,false);
        apply_l2p_flags( l2p_mask, bt08, j,false);
        apply_l2p_flags( l2p_mask, sst, j,false);
        time_ind++;
    }


    filename = generate_filename(paths[file_count]);
    generate_masks(bt08,bt10,bt11,bt12, sst,reference, border_mask, bt_mask,current_mid);
    

    filename = generate_filename(paths[file_count]);
    file_count++;

    while(time_ind < sample_size+1){


        save_loc = folder_loc+ filename;       
        save_test_nc_float(bt_mask,save_loc.c_str());

        clear_paths.push_back(save_loc.c_str());
        printf("Generated File %s\n",save_loc.c_str());
        
        if(time_ind < sample_size){
            printf("continuing loop\n");
            readgranule(paths[time_ind], bt11,bt12,bt08,bt10, sst, new_granule_ind);
            printf("read new granule\n");

            
            apply_l2p_flags( l2p_mask, bt11, new_granule_ind,false);
            apply_l2p_flags( l2p_mask, bt12, new_granule_ind,false);
            apply_l2p_flags( l2p_mask, bt10, new_granule_ind,false);
            apply_l2p_flags( l2p_mask, bt08, new_granule_ind,false);
            apply_l2p_flags( l2p_mask, sst, new_granule_ind,false);
            new_granule_ind = (new_granule_ind+1) % FILTER_TIME_SIZE;
            printf("starting mask generation\n");
            
            filename = generate_filename(paths[file_count]);
            generate_masks(bt08,bt10,bt11,bt12,sst, reference, border_mask, bt_mask,new_mask_ind);
            new_mask_ind = (new_mask_ind+1) % FILTER_TIME_SIZE;  
            current_mid =(current_mid + 1) % FILTER_TIME_SIZE; 
            
            file_count++;         
        }
        time_ind++;i++;
    }
    
}
       
    