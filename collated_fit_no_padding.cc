void
calc_approximate(const Mat1f &bt11_smooth, const Mat1f &bt11_clear, const Mat1b &l2p_mask, Mat1f &bt11_approx, int time_size)
{
    int x,y,z;
    double d1,d2,coeff;
    //placeholder
    int count=0;
    vector<float> data_cl2;
    vector<float> smooth;
    vector<int> clear_inds;
    //Mat1f diffs(HEIGHT,WIDTH);
    bt11_approx.setTo(NAN);
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < WIDTH; ++x){
            if(l2p_mask(y,x) == 0 ){
                count = 0;
                d1=0;d2 =0;
                data_cl2.clear();
                smooth.clear();
                clear_inds.clear();
                for(z = 0; z < time_size; ++z){
                    
                    if(!std::isnan(bt11_smooth(y,x,z)) && !std::isnan(bt11_clear(y,x,z))){
                        //if(!std::isnan(bt11_clear(y,x,z)) || (fabs(bt11_smooth(y,x,z) - sst(y,x,z))< T_SST_SMOOTH)){

                        d1 += bt11_clear(y,x,z)*bt11_smooth(y,x,z);
                        d2 += bt11_smooth(y,x,z)*bt11_smooth(y,x,z);
                        count++;
                                
                    }                   
                }
            
                //diffs(y,x) = d2-d1;
                if(d2 != 0 && count > MIN_POINTS){
                    coeff = d1/d2;
                    for(z = 0; z < time_size; ++z){
                        bt11_approx(y,x,z)= coeff * bt11_smooth(y,x,z);
                    }
                }
            }
        }
    }
}

void
collate_samples_smooth(const Mat1f &approx, Mat1f &collated, const MatrixXf &Gamma,  
                const vector<int> collated_inds, const int sample_size, const Mat1b &l2p_mask)
{
    int y,x,i,j,k,t,ind;
    int collated_size = collated_inds.size();
    int w = 1;
    int first_valid = -1;
    int last_valid = -1;
    float approx_count, approx_sum;
    VectorXf collated_vals(collated_size);
    VectorXf rhs(collated_size);
    MatrixXf v(collated_size,collated_size);

    rhs.setZero(); v.setZero();

    //for each pixel
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < HEIGHT; ++x){
            first_valid = -1; last_valid = -1;
            if(l2p_mask(y,x) ==0){
                //for each hour in collated
                //auto start = std::chrono::system_clock::now();
                for(i = 0; i < collated_size; ++i){
                    //for each value in current collated hour
                    approx_sum = approx_count = 0;
                    ind = collated_inds[i];
                    
                    // handle case where there are not enough granules in the front of the interval 
                    if(ind < COLLATED_LAG){
                        for(t = -ind;t < COLLATED_LAG+1; ++t){
                            for(k = -w; k < w+1; ++k){
                                for(j = -w; j < w+1; ++j){                                
                                    if(std::isfinite(approx(y+k,x+j,ind+t))){
                                        approx_sum += approx(y+k,x+j,ind+t);
                                        approx_count += 1;
                                    }
                                }
                            }                            
                        }
                    }

                    else if(ind >= COLLATED_LAG && ((sample_size - ind) > COLLATED_LAG)){
                        for(t = -COLLATED_LAG; t < COLLATED_LAG+1; ++t){
                            for(k = -w; k < w+1; ++k){
                                for(j = -w; j < w+1; ++j){                                    
                                    if(std::isfinite(approx(y+k,x+j,ind+t))){
                                        approx_sum += approx(y+k,x+j,ind+t);
                                        approx_count += 1;
                                    }
                                }
                            }                            
                        }
                    }

                    else if(sample_size - ind <= COLLATED_LAG){
                        for(t = -COLLATED_LAG; t < (sample_size - ind); ++t){
                            for(k = -w; k < w+1; ++k){
                                for(j = -w; j < w+1; ++j){                                    
                                    if(std::isfinite(approx(y+k,x+j,ind+t))){
                                        approx_sum += approx(y+k,x+j,ind+t);
                                        approx_count += 1;
                                    }
                                }
                            }                            
                        }
                    }

                    else eprintf("collate_smooth: index %d invalid\n", ind);
                    

                    /*
                    if(approx_count < win_size/2){
                        approx_count = 0;
                        approx_sum = 0;
                    }
                    */
                    rhs[i] = MU_APPROX*approx_sum ;
                    v(i,i) = MU_APPROX*approx_count;
                    if(v(i,i)!=0){                      
                        if(first_valid == -1){
                            first_valid = i;
                        }
                        if(last_valid < i){
                            last_valid = i;
                        }
                    }
                }
                

                if(last_valid != -1 && first_valid != -1){
                    collated_vals = (v + Gamma).llt().solve(rhs);              
                    for(i=first_valid;i<=last_valid;i++){
                        collated(y,x,i) = collated_vals(i);
                    }                  
                }
            }
        }
    }
}

void
collate_samples(const Mat1f &clear_samples, const Mat1f &approx, Mat1f &collated, const MatrixXf &Gamma,  
    const vector<int> collated_inds,const Mat1b &l2p_mask, int sample_size, const string ref_file)
{
    
    // 1) set up right hand side-> rhs[k] = MU*clear_sum + MU* approx sum -> k = time step
    // 2) set up left hand side -> diagonal = mu*clear_count+mu*approx_count
    // 3) add gamma to left hand side for tridiagonal matrix
    // 4) these should all be premade eigen arrays determined by the number of collated values being produced
    // 5) collated = inv(lhs)*rhs

    // read reference for sst_reynolds and o2ld
    // if abs(collated - reynolds ) > 7 && o2ld <= 50 -> sset to nan
    int y,x,i,t,ind;
    int collated_size = collated_inds.size();
    int half_ind = COLLATED_LAG;
    int first_valid = -1;
    int last_valid = -1;
    float clear_sum, clear_count, approx_count, approx_sum;
    VectorXf collated_vals(collated_size);
    VectorXf rhs(collated_size);
    MatrixXf v(collated_size,collated_size);

    Mat1f reference(HEIGHT,WIDTH);
    Mat1f o2ld(HEIGHT,WIDTH);

    rhs.setZero(); v.setZero();

    get_var(ref_file,reference,"sst_reynolds");
    get_var(ref_file,o2ld,"ocean_to_land_dist");

    //for each pixel
    for(y = 0; y < HEIGHT; ++y){
        for(x = 0; x < HEIGHT; ++x){
            first_valid = -1; last_valid = -1;
            if(l2p_mask(y,x) ==0){
                //for each hour in collated
                //auto start = std::chrono::system_clock::now();
                for(i = 0; i < collated_size; ++i){
                    //for each value in current collated hour
                    clear_sum = clear_count = approx_sum = approx_count = 0;
                    ind = collated_inds[i]; 
                    
                    if(ind < COLLATED_LAG){
                        for(t = -ind; t < COLLATED_LAG+1; ++t){
                            if(std::isfinite(clear_samples(y,x,ind+t))){
                                clear_sum += clear_samples(y,x,ind+t);
                                clear_count++;
                            }
                            if(std::isfinite(approx(y,x,ind+t))){
                                approx_sum += approx(y,x,ind+t);
                                approx_count++;
                            }                       
                        }
                        half_ind = (ind + COLLATED_LAG) /2;

                    }

                    else if(ind >= COLLATED_LAG && ((sample_size - ind) > COLLATED_LAG)){
                        for(t = -COLLATED_LAG; t < COLLATED_LAG+1; ++t){
                            if(std::isfinite(clear_samples(y,x,ind+t))){
                                clear_sum += clear_samples(y,x,ind+t);
                                clear_count++;
                            }
                            if(std::isfinite(approx(y,x,ind+t))){
                                approx_sum += approx(y,x,ind+t);
                                approx_count++;
                            }                    
                        }
                        half_ind = COLLATED_LAG;
                    }

                    else if(sample_size - ind <= COLLATED_LAG){
                        for(t = -COLLATED_LAG; t < (sample_size - ind); ++t){
                            if(std::isfinite(clear_samples(y,x,ind+t))){
                                clear_sum += clear_samples(y,x,ind+t);
                                clear_count++;
                            }
                            if(std::isfinite(approx(y,x,ind+t))){
                                approx_sum += approx(y,x,ind+t);
                                approx_count++;
                            }                    
                            half_ind = (COLLATED_LAG + sample_size - ind)/2;
                        }
                    }

                    else eprintf("collate_smooth: index %d invalid\n", ind);

                    if(clear_count < 3 && approx_count == 0){
                        clear_count =0;
                        clear_sum = 0;
                    }
                    if(approx_count < half_ind){
                        approx_count = 0;
                        approx_sum = 0;
                    }
                    
                    rhs[i] = MU_APPROX*approx_sum + MU_CLEAR*clear_sum;
                    v(i,i) = MU_APPROX*approx_count +MU_CLEAR*clear_count;
                    if(v(i,i)!=0){                      
                        if(first_valid == -1){
                            first_valid = i;
                        }
                        if(last_valid < i){
                            last_valid = i;
                        }
                    }
                }
                

                if(last_valid != -1 && first_valid != -1){

                    collated_vals = (v + Gamma).llt().solve(rhs);
              
                    for(i = first_valid;i <= last_valid; ++i){
                        collated(y,x,i) = collated_vals(i);
                        if(abs(reference(y,x) - collated(y,x,i)) > 7 && o2ld(y,x) <=50){
                            collated(y,x,i) = NAN;
                        }
                    }
                          
                }
            }
        }
    }

}

void
approx_clear(const vector<string> &smooth_paths, const vector<string> &clear_paths, const vector<string> &original_paths,
             const Mat1b &l2p_mask, const string ref_file, bool interp)
{    


    int sample_size = smooth_paths.size();
    //int sample_size = 50;
    //interpolate checks the n+1 index to ensure their derivative is small
    //approx is HEIGHT x WIDTH x time_size-1
    //int time_size = 10;
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG;

    int i,j,y,x;
    int dims[3] = {HEIGHT,WIDTH,sample_size};


    Mat1f smooth_samples(3,dims);
    Mat1f approx_samples(3,dims);
    approx_samples.setTo(NAN);

    string approx1_folder_loc = "data/approx";
    string collated1_folder_loc = "data/collated_mat";
    string approx2_folder_loc = "data/approx2";
    string collated2_folder_loc = "data/collated_mat2";
    string smooth_folder_loc = "data/smooth_collate";
    string reinstated_folder_loc = "data/reinstated";
    string filename, save_loc;
    
    vector<int> collated_inds; //indicies where collation takes place
    vector<int> interp_inds; //indicies where collated values sit in interpolation matrix

    // save locations for generated data
    vector<string> smooth2_paths; // collated smooth
    vector<string> approx1_paths;
    vector<string> collated1_paths;
    vector<string> approx2_paths;
    vector<string> collated2_paths;
    vector<string> reinstated_paths;
    vector<string> collated_hour_paths;

    /////////////////////////////////////////////////
    // FIRST COMPUTE FILE NAMES AND SIZES/INDICIES //
    // OF OUTPUT MATRICES                          //
    /////////////////////////////////////////////////

    //generate save paths for collated_values as well as indicies of hour
    for(i = 0; i < sample_size - 1; ++i){
        filename = generate_filename(smooth_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            collated_inds.push_back(i);
            collated_hour_paths.push_back(collated2_folder_loc + filename);
        }
    }

    int collated_size = collated_inds.size(); // how many hours of colaltion we have
    printf("collated size = %d\n",collated_size);
    // filenames and indicies for the rest of the interval
    // last granule is ignored to account for the removal of large derivatives
    for(i = collated_inds[0]  ; i <= collated_inds[collated_size - 1]; ++i){
        filename = generate_filename(smooth_paths[i]);
    
        smooth2_paths.push_back(smooth_folder_loc + filename);
        approx1_paths.push_back(approx1_folder_loc + filename);
        collated1_paths.push_back(collated1_folder_loc + filename);
        approx2_paths.push_back(approx2_folder_loc + filename);
        collated2_paths.push_back(collated2_folder_loc + filename);
        reinstated_paths.push_back(reinstated_folder_loc + filename);
    }
    
    // Init matricies to perform collation steps
   
    int collated_interp_size = smooth2_paths.size(); // size after interpolating hourly collation
    int dims_collated[3] ={HEIGHT,WIDTH,collated_size};
    int dims_interpolated[3] = {HEIGHT,WIDTH,collated_interp_size};

    
    Mat1f collated(3,dims_collated);
    Mat1f collated_interp(3,dims_interpolated);
    collated.setTo(NAN);
    collated_interp.setTo(NAN);

    // calculate the x positions for the interpolation steps, all x's are on the hour
    // however since granules may be missing, distance between hours may be variable
    for(i = 0; i < collated_interp_size; ++i){
        filename = generate_filename(smooth2_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            interp_inds.push_back(i);
        }
    }

    
    // Init Gamma matrix for collation procedure
    MatrixXf Gamma(collated_size,collated_size);
    Gamma.setZero();

    //fill gamma matrix
    //first and last row have different values
    Gamma(0,0) = GAMMA;Gamma(0,1) = -GAMMA;
    Gamma(collated_size-1,collated_size-1) = GAMMA; Gamma(collated_size-1,collated_size-2) = -GAMMA;
    for(i=1;i<collated_size-1;i++){
        Gamma(i,i)   = 2*GAMMA;
        Gamma(i,i-1) = -GAMMA;
        Gamma(i,i+1) = -GAMMA;
    }
    
    printf("GAMMA initializaed\n");
    // open all smooth files
    printf("getting ice masks and cloud masks\n");
    for(j=0;j<sample_size;j++){
        readgranule_oneband(smooth_paths[j],smooth_samples,j,"sea_surface_temperature");
        printf("reading file %d for variable %s\n",j,"sea_surface_temperature");
    }

    // remove pixels with high derivatives
    printf("starting removing derivatives\n");
    remove_high_derivatives(smooth_samples, l2p_mask, sample_size);
    sample_size--; // handle removal of final slice
    printf("Finished removing derivatives\n");
    
    auto start = std::chrono::system_clock::now();
    // collate smooth samples
    // parameter approx samples is all NAN and will not count towards collation
    printf("Starting Collation on sea_surface_temperature\n");
    collate_samples_smooth(smooth_samples, collated, Gamma, collated_inds, sample_size, l2p_mask);
    printf("Finished Collation on sea_surface_temperature\n");
    
    smooth_samples.release();
    //save_mat(smooth2_paths, smooth_collated, "sea_surface_temperature",true);

    for(i = 0 ; i < collated_size; ++i){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_interp(y,x,interp_inds[i]) = collated(y,x,i);
            }
        }
    }
    interpolate_hermite(collated_interp,l2p_mask,collated_interp_size,T_INTERP, 2*INTERP_DIST,false);
    auto end = std::chrono::system_clock::now();
    auto elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to open clear files = " << elapsed.count() << '\n';
    
    ////////////////////////////////////////
    // APPROXIMATION WITH COLLATED SMOOTH //
    // AND FIRST COLLATING PASS           //
    ////////////////////////////////////////
    Mat1b clear_masks(HEIGHT,WIDTH);
    Mat1f clear_samples(3,dims_interpolated);
    Mat1f approx(3,dims_interpolated);

    start = std::chrono::system_clock::now();
    for(j = 0; j < collated_interp_size; ++j){
        read_mask(clear_paths[j+collated_inds[0]],clear_masks,-1);
        readgranule_oneband(original_paths[j+original_lag+collated_inds[0]],clear_samples,j,"sea_surface_temperature");
        apply_mask_slice(clear_masks,clear_samples,j,false);
        printf("read file1 %s\n",original_paths[j+original_lag+collated_inds[0]].c_str());
    }
    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to open clear files = " << elapsed.count() << '\n';

    start = std::chrono::system_clock::now();
    printf("calculating interp clear\n");

    interpolate_hermite(clear_samples, l2p_mask, collated_interp_size, T_INTERP, 2*INTERP_DIST, false);

    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to open interpolate clear = " << elapsed.count() << '\n';

    start = std::chrono::system_clock::now();
    printf("calculating approximation for sea_surface_temperature\n");
    calc_approximate(collated_interp,clear_samples, l2p_mask, approx, collated_interp_size);
    printf("finished approximation for sea_surface_temperature\n");
    //save_mat(approx1_paths, approx, "sea_surface_temperature",true);
    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to open finish approx 1 = " << elapsed.count() << '\n';

    start = std::chrono::system_clock::now();
    // interp inds are now the location of the hour in for the clear samples and approx matrices
    printf("Starting Collation on sea_surface_temperature\n");
    collated.setTo(NAN);
    collate_samples(clear_samples, approx, collated, Gamma,interp_inds,l2p_mask, collated_interp_size, ref_file);
    printf("Finished Collation on sea_surface_temperature\n");
    
    //SAVENC(collated);
    approx.release();
    //clear_samples.release();
    
    collated_interp.setTo(NAN);
    for(i = 0; i < collated_size; ++i){
        for(y = 0; y < HEIGHT; ++y){
            for( x = 0; x < WIDTH; ++x){
                collated_interp(y,x,interp_inds[i]) = collated(y,x,i);
            }
        }        
    }

    printf("starting interpolation of collated values\n");
    interpolate_hermite(collated_interp,l2p_mask,collated_interp_size,T_INTERP,2*INTERP_DIST,false);
    //remove_last_value(collated_interp,invalid_mask,land_mask,collated_interp_size);

    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to finish first collation = " << elapsed.count() << '\n';

    printf("finished collated interpolation\n");
    //save_mat(collated1_paths, collated_interp, "sea_surface_temperature",true);
    
    /////////////////////////////////////////
    // SMOOTH INTERPOLATED COLLATED VALUES //
    /////////////////////////////////////////
    
    int collated_smooth_lag = 3;    
    int window[3] = {5,5,collated_smooth_lag};

    Mat1f collated_smooth(3,dims_interpolated);
    Mat1f reinstated_clear(3,dims_interpolated);
    Mat1f original_sst(HEIGHT,WIDTH);
    float DD;
    string full_path;

    start = std::chrono::system_clock::now();
    printf("starting smoothing of collated values\n");
    smooth_samples_collated(collated_interp, collated_smooth,l2p_mask, ref_file, window, collated_interp_size);
    printf("finished smoothing of collated values\n");
    //save_mat(smooth2_paths, collated_smooth, "sea_surface_temperature",true);
    
    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to finish smooth collated = " << elapsed.count() << '\n';

    /////////////////////////////////////////////////
    // REINSTATE CLEAR WHERE ORIGINAL SST > SMOOTH //
    /////////////////////////////////////////////////

    //get folder where original files are located
    Mat1f clear_sst(HEIGHT,WIDTH);

    start = std::chrono::system_clock::now();
    for(i = 0; i < collated_interp_size; ++i){        
        get_var(original_paths[i+original_lag+collated_inds[0]], original_sst, "sea_surface_temperature");
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                clear_sst(y,x) = original_sst(y,x);
            }
        }
        read_mask(clear_paths[i+collated_inds[0]],clear_masks,-1);
        apply_mask_2d(clear_masks, clear_sst);

        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                reinstated_clear(y,x,i) = clear_sst(y,x);
                DD = collated_smooth(y,x,i) - original_sst(y,x);
                if(std::isfinite(DD) && DD < T_SMOOTH_COLLATED){
                    reinstated_clear(y,x,i) = original_sst(y,x);
                }
            }
        }
    }
    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to finish smooth collated = " << elapsed.count() << '\n';

    clear_sst.release();
    original_sst.release();
    clear_masks.release();
    save_mat(reinstated_paths, reinstated_clear, "sea_surface_temperature",true);

    ///////////////////////////////////////////////
    // APPROXIMATE AND COLLATE FINAL TIME        //
    // WITH REINSTATED CLEAR AND SMOOTH COLLATED //
    ///////////////////////////////////////////////
    
    clear_samples.release();
    //approximate using new clear and smooth collated
    //collate one last time
    
    start = std::chrono::system_clock::now();
    approx.create(3,dims_interpolated);
    printf("calculating approximation for sea_surface_temperature\n");
    calc_approximate(collated_smooth,reinstated_clear, l2p_mask, approx, collated_interp_size);
    end = std::chrono::system_clock::now();
    elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "time to finish approx2 = " << elapsed.count() << '\n';

    printf("finished approximation for sea_surface_temperature\n");
    //Mat1f collated_approx(3,dims);
    save_mat(approx2_paths, approx, "sea_surface_temperature",true);


    //start = std::chrono::system_clock::now();
    collated.setTo(NAN);
    printf("Starting Collation on sea_surface_temperature\n");
    collate_samples(reinstated_clear, approx, collated, Gamma, interp_inds,l2p_mask,collated_interp_size, ref_file);
    printf("Finished Collation on sea_surface_temperature\n");
    reinstated_clear.release();approx.release();

    


    collated_interp.setTo(NAN);
    for(i = 0; i < collated_size; ++i){
        for(y = 0; y < HEIGHT; ++y){
            for(x = 0; x < WIDTH; ++x){
                collated_interp(y,x,interp_inds[i]) = collated(y,x,i);
            }
        }
        
    }
    
    //collated.release();
    printf("starting interpolation of collated values\n");
    if(interp){
        interpolate_hermite(collated_interp,l2p_mask,collated_interp_size,T_INTERP,2*INTERP_DIST,false);
        remove_last_value(collated_interp,collated_interp_size);
        save_mat(collated2_paths, collated_interp, "sea_surface_temperature",true);
    }
    //end = std::chrono::system_clock::now();
    //elapsed =  std::chrono::duration_cast<std::chrono::seconds>(end - start);
    //std::cout << "time to finish second collation = " << elapsed.count() << '\n';
    else{
        save_mat(collated_hour_paths, collated, "sea_surface_temperature",true);
    }
    collated.release();
    collated_interp.release();
    
}

