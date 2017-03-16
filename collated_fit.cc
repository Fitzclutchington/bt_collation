void
calc_approximate(const Mat1f &bt11_smooth, const Mat1f &bt11_clear, const Mat1b &land_mask, const Mat1b &invalid_mask, 
                 Mat1f &bt11_approx, int time_size)
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
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(land_mask(y,x) == 0 && invalid_mask(y,x) == 0 ){
                count = 0;
                d1=0;d2 =0;
                data_cl2.clear();
                smooth.clear();
                clear_inds.clear();
                for(z=0;z<time_size;z++){
                    
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
                    for(z=0;z<time_size;z++){
                        bt11_approx(y,x,z)= coeff * bt11_smooth(y,x,z);
                    }
                }
            }
        }
    }
}

void
collate_samples(const Mat1f &approx, Mat1f &collated, const MatrixXf &Gamma,  
    const vector<int> collated_inds,const Mat1b &land_mask,const Mat1b &invalid_mask, bool interp)
{
    //TODO:
    // 1) set up right hand side-> rhs[k] = MU*clear_sum + MU* approx sum -> k = time step
    // 2) set up left hand side -> diagonal = mu*clear_count+mu*approx_count
    // 3) add gamma to left hand side for tridiagonal matrix
    // 4) these should all be premade eigen arrays determined by the number of collated values being produced
    // 5) collated = inv(lhs)*rhs
    int y,x,i,j,k,t,ind;
    int collated_size = collated_inds.size();
    int w = 1;
    int win_size = (2*COLLATED_LAG +1)*(2*w +1)*(2*w+1);
    int first_valid = -1;
    int last_valid = -1;
    float clear_sum, clear_count, approx_count, approx_sum;
    VectorXf collated_vals(collated_size);
    VectorXf rhs(collated_size);
    MatrixXf v(collated_size,collated_size);

    rhs.setZero(); v.setZero();

    //for each pixel
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<HEIGHT;x++){
            first_valid = -1; last_valid = -1;
            if(land_mask(y,x)==0 && invalid_mask(y,x) ==0){
                //for each hour in collated
                //auto start = std::chrono::system_clock::now();
                for(i=0;i<collated_size;i++){
                    //for each value in current collated hour
                    clear_sum = clear_count = approx_sum = approx_count = 0;
                    ind = collated_inds[i]; 
                    //TODO: add window 3x3
                    for(t=-COLLATED_LAG;t<COLLATED_LAG+1;t++){
                        for(k=-w;k<w+1;k++){
                            for(j=-w;j<w+1;j++){
                                
                                if(!std::isnan(approx(y+k,x+j,ind+t))){
                                    approx_sum += approx(y+k,x+j,ind+t);
                                    approx_count += 1;
                                }
                            }
                        }
                        
                    }

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
                    //auto end = std::chrono::system_clock::now();
                    //auto elapsed =  std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    //std::cout << "time to compute rhs = " << elapsed.count() << '\n';

                    //start = std::chrono::system_clock::now();
                    collated_vals = (v + Gamma).llt().solve(rhs);
                    //end = std::chrono::system_clock::now();
                    //elapsed =  std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    //std::cout << "time to compute multiplication = " << elapsed.count() << '\n';              
                    for(i=first_valid;i<=last_valid;i++){
                        collated(y,x,i) = collated_vals(i);
                    }
                    //collated(y,x,first_valid) = NAN; collated(y,x,last_valid)= NAN;                   
                }
            }
        }
    }

}

void
collate_samples_final(const Mat1f &clear_samples, const Mat1f &approx, Mat1f &collated, const MatrixXf &Gamma,  
    const vector<int> collated_inds,const Mat1b &land_mask,const Mat1b &invalid_mask, bool interp, const string ref_file)
{
    //TODO:
    // 1) set up right hand side-> rhs[k] = MU*clear_sum + MU* approx sum -> k = time step
    // 2) set up left hand side -> diagonal = mu*clear_count+mu*approx_count
    // 3) add gamma to left hand side for tridiagonal matrix
    // 4) these should all be premade eigen arrays determined by the number of collated values being produced
    // 5) collated = inv(lhs)*rhs

    // read reference for sst_reynolds and o2ld
    // if abs(collated - reynolds ) > 7 && o2ld <= 50 -> sset to nan
    int y,x,i,j,k,t,ind;
    int collated_size = collated_inds.size();
    int w = 1;
    int win_size = (2*COLLATED_LAG +1)*(2*w +1)*(2*w+1);
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
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<HEIGHT;x++){
            first_valid = -1; last_valid = -1;
            if(land_mask(y,x)==0 && invalid_mask(y,x) ==0){
                //for each hour in collated
                //auto start = std::chrono::system_clock::now();
                for(i=0;i<collated_size;i++){
                    //for each value in current collated hour
                    clear_sum = clear_count = approx_sum = approx_count = 0;
                    ind = collated_inds[i]; 
                    //TODO: add window 3x3
                    for(t=-COLLATED_LAG;t<COLLATED_LAG+1;t++){
                        if(!std::isnan(clear_samples(y,x,ind+t))){
                            clear_sum += clear_samples(y,x,ind+t);
                            clear_count++;
                        }
                        if(!std::isnan(approx(y,x,ind+t))){
                            approx_sum += approx(y,x,ind+t);
                            approx_count += 1;
                        }                    
                        
                    }

                    if(clear_count < 3 && approx_count == 0){
                        clear_count =0;
                        clear_sum = 0;
                    }
                    if(approx_count < COLLATED_LAG){
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
                    //auto end = std::chrono::system_clock::now();
                    //auto elapsed =  std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    //std::cout << "time to compute rhs = " << elapsed.count() << '\n';

                    //start = std::chrono::system_clock::now();
                    collated_vals = (v + Gamma).llt().solve(rhs);
                    //end = std::chrono::system_clock::now();
                    //elapsed =  std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                    //std::cout << "time to compute multiplication = " << elapsed.count() << '\n';              
                    for(i=first_valid;i<=last_valid;i++){
                        collated(y,x,i) = collated_vals(i);
                        if(abs(reference(y,x) - collated(y,x,i)) > 7 && o2ld(y,x) <=50){
                            collated(y,x,i) = NAN;
                        }
                    }
                    //collated(y,x,first_valid) = NAN; collated(y,x,last_valid)= NAN;                   
                }
            }
        }
    }

}

void
approx_clear(const vector<string> &smooth_paths, const vector<string> &clear_paths, const vector<string> &original_paths,
             const Mat1b &land_mask, const Mat1b &invalid_mask, vector<string> &approx_paths, const string ref_file, bool interp){    


    int time_size = smooth_paths.size();
    //int time_size = 60;
    //interpolate checks the n+1 index to ensure their derivative is small
    //approx is HEIGHT x WIDTH x time_size-1
    int approx_time = time_size-1;
    //int time_size = 10;
    int original_lag = FILTER_WINDOW_LAG+SECOND_PASS_LAG+SMOOTH_WINDOW_LAG;

    int i,j,y,x;
    int dims[3] = {HEIGHT,WIDTH,time_size};
    int dims_approx[3] = {HEIGHT,WIDTH,time_size-1};


    Mat1f smooth_samples(3,dims);

    string approx_folder_loc = "data/approx";
    string collated_folder_loc = "data/collated_mat";
    string smooth_folder_loc = "data/smooth_collate";
    string filename, save_loc;
    
    vector<int> collated_inds;
    vector<int> interp_inds;
    vector<int> smooth_inds;

    vector<string> collated_paths;
    vector<string> collated_interp_paths;
    vector<string> smooth_collate;
    vector<string> smooth2_paths;


    //generate save paths for collated_values as well as indicies of hour
    for(i=COLLATED_LAG;i<time_size - COLLATED_LAG;i++){
        filename = generate_filename(smooth_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            smooth_collate.push_back(filename);
            smooth_inds.push_back(i);
        }
    }

    
    //generate filenames for interp collated
    for(i=smooth_inds[0];i<time_size - COLLATED_LAG;i++){
        filename = generate_filename(smooth_paths[i]);  
        save_loc = smooth_folder_loc + filename;
        smooth2_paths.push_back(save_loc);
    }


    
    int smooth_size = smooth_collate.size();
    int smooth_interp_size = smooth2_paths.size();
    int dims_smooth2[3] ={HEIGHT,WIDTH,smooth_size};

    for(i=0;i<smooth_interp_size;i++){
        filename = generate_filename(smooth2_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            interp_inds.push_back(i);
        }
    }
    
    Mat1f smooth_collated(3,dims_smooth2);
    smooth_collated.setTo(NAN);
    MatrixXf Gamma(smooth_size,smooth_size);
    Gamma.setZero();

    //fill gamma matrix
    //first and last row have different values
    Gamma(0,0) = GAMMA;Gamma(0,1) = -GAMMA;
    Gamma(smooth_size-1,smooth_size-1) = GAMMA; Gamma(smooth_size-1,smooth_size-2) = -GAMMA;
    for(i=1;i<smooth_size-1;i++){
        Gamma(i,i)   = 2*GAMMA;
        Gamma(i,i-1) = -GAMMA;
        Gamma(i,i+1) = -GAMMA;
    }

    printf("getting ice masks and cloud masks\n");
    for(j=0;j<time_size;j++){
        //read_mask(clear_paths[j+SMOOTH_WINDOW_LAG],clear_masks,j);
        //read_acspo(original_paths[j+original_lag],clear_masks,j);
        //readgranule_oneband(original_paths[j+original_lag],clear_samples,j,"sea_surface_temperature");
        //apply_mask_slice(clear_masks,clear_samples,j,true);
        readgranule_oneband(smooth_paths[j],smooth_samples,j,"sea_surface_temperature");
        printf("reading file %d for variable %s\n",j,"sea_surface_temperature");
    }

    //TODO WRITE THIS FUNCTION 
    printf("starting removing derivatives\n");
    remove_high_derivatives(smooth_samples, land_mask, invalid_mask, time_size);
    printf("Finished removing derivatives\n");

    printf("Starting Collation on brightness_temperature_08um6\n");
    collate_samples(smooth_samples, smooth_collated, Gamma,  smooth_inds,land_mask,invalid_mask,interp);
    printf("Finished Collation on brightness_temperature_08um6\n");
    
    smooth_samples.release();
    //save_mat(smooth2_paths, smooth_collated, "sea_surface_temperature",true);

    dims[2] = smooth_interp_size;
    Mat1f smooth_interp(3,dims);
    smooth_interp.setTo(NAN);
    for(i=0;i<interp_inds.size();i++){
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                smooth_interp(y,x,interp_inds[i]) = smooth_collated(y,x,i);
            }
        }
        
    }
    interpolate_hermite(smooth_interp,land_mask,invalid_mask,smooth_interp_size,T_INTERP, 2*INTERP_DIST,false);
    //save_mat(smooth2_paths, smooth_interp, "sea_surface_temperature",true);
    
    smooth_collated.release();

    Mat1b clear_masks(HEIGHT,WIDTH);
    Mat1f clear_samples(3,dims);
    Mat1f approx(3,dims);

    printf("smooth inds lag = %d\n",smooth_inds[0]);

    for(j=0;j<smooth_interp_size;j++){
        read_mask(clear_paths[j+SMOOTH_WINDOW_LAG+smooth_inds[0]],clear_masks,-1);
        readgranule_oneband(original_paths[j+original_lag+smooth_inds[0]],clear_samples,j,"sea_surface_temperature");
        apply_mask_slice(clear_masks,clear_samples,j,false);
        printf("read file %s\n",original_paths[j+original_lag+smooth_inds[0]].c_str());
    }

    printf("starting generation of approx filenames\n");
    for(j=0;j<smooth_interp_size;j++){
        filename = generate_filename(smooth2_paths[j]);
        save_loc = approx_folder_loc+filename;
        printf("filename = %s\n", save_loc.c_str());
        approx_paths.push_back(save_loc);
    }

    printf("starting generation of collated\n");
    //generate save paths for collated_values as well as indicies of hour
    // previously < smooth_interp_size
    for(i=COLLATED_LAG;i<approx_paths.size() - COLLATED_LAG;i++){
        filename = generate_filename(approx_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            collated_paths.push_back(filename);
            collated_inds.push_back(i);
        }
    }
    
    printf("starting generation of collated interp\n");
    //generate filenames for interp collated
    for(i=collated_inds[0];i<approx_paths.size()-COLLATED_LAG;i++){
        filename = generate_filename(approx_paths[i]); 
        save_loc = collated_folder_loc + filename; 
        printf("approx filename = %s\n",save_loc.c_str());    
        collated_interp_paths.push_back(save_loc);
    }

    int collated_size = collated_paths.size();
    int collated_interp_size = collated_interp_paths.size();
    int dims_collated[3] ={HEIGHT,WIDTH,collated_size};

    printf("starting collated interp \n");
    interp_inds.clear();
    for(i=0;i<collated_interp_size;i++){
        printf("filename = %s\n",collated_interp_paths[i].c_str());
        filename = generate_filename(collated_interp_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            interp_inds.push_back(i);
        }
    }
    
    printf("starting gamma init\n"); 
    Mat1f collated(3,dims_collated);
    collated.setTo(NAN);
    Gamma.resize(collated_size,collated_size);
    Gamma.setZero();

    printf("(starting GAMMA\n" );
    //fill gamma matrix
    //first and last row have different values
    Gamma(0,0) = GAMMA;Gamma(0,1) = -GAMMA;
    Gamma(collated_size-1,collated_size-1) = GAMMA; Gamma(collated_size-1,collated_size-2) = -GAMMA;
    for(i=1;i<collated_size-1;i++){
        Gamma(i,i)   = 2*GAMMA;
        Gamma(i,i-1) = -GAMMA;
        Gamma(i,i+1) = -GAMMA;
    }

    printf("calculating interp clear\n");
    interpolate_hermite(clear_samples, land_mask, invalid_mask, smooth_interp_size, T_INTERP, 2*INTERP_DIST, false);

    printf("calculating approximation for brightness_temperature_08um6\n");
    calc_approximate(smooth_interp,clear_samples, land_mask,invalid_mask, approx,smooth_interp_size);
    printf("finished approximation for brightness_temperature_08um6\n");
    smooth_interp.release();

    printf("Starting Collation on brightness_temperature_08um6\n");
    collate_samples_final(clear_samples, approx, collated, Gamma,  collated_inds,land_mask,invalid_mask,interp,ref_file);
    printf("Finished Collation on brightness_temperature_08um6\n");
    //save_mat(approx_paths, approx, "sea_surface_temperature",true);
    approx.release();
    //clear_samples.release();
    
    dims[2] = collated_interp_size;
    Mat1f collated_interp(3,dims);
    collated_interp.setTo(NAN);
    for(i=0;i<interp_inds.size();i++){
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                collated_interp(y,x,interp_inds[i]) = collated(y,x,i);
            }
        }
        
    }
    collated.release();
    printf("starting interpolation of collated values\n");
    interpolate_hermite(collated_interp,land_mask,invalid_mask,collated_interp_size,T_INTERP,2*INTERP_DIST,false);
    //remove_last_value(collated_interp,invalid_mask,land_mask,collated_interp_size);

    printf("finished collated interpolation\n");
    //save_mat(collated_interp_paths, collated_interp, "sea_surface_temperature",true);


    int collated_smooth_lag = 3;    
    int window[3] = {5,5,collated_smooth_lag};
    dims[0] = HEIGHT;dims[1] = WIDTH; dims[2] = collated_interp_paths.size()-6;
    Mat1f collated_smooth(3,dims);
    Mat1f reinstated_clear(3,dims);
    Mat1f original_sst(HEIGHT,WIDTH);
    float DD;
    string full_path;
    vector<string> collated_smooth_paths;
    vector<string> collated_approx_paths;

    for(i=3;i<collated_interp_paths.size()-3;i++){
        filename = generate_filename(collated_interp_paths[i]);
        save_loc = approx_folder_loc + filename;
        collated_approx_paths.push_back(save_loc);
        collated_smooth_paths.push_back(collated_interp_paths[i]);
    }
    smooth_samples_collated(collated_interp, collated_smooth,land_mask, invalid_mask, ref_file, window, collated_interp_paths.size());
    
    //get folder where original files are located
    string original_folder_loc = generate_foldername(original_paths[0]);
    for(i=0;i<collated_smooth_paths.size();i++){
        original_folder_loc = generate_foldername(original_paths[i+original_lag+smooth_inds[0]+collated_inds[0]+collated_smooth_lag]);
        filename = generate_filename(collated_smooth_paths[i]);
        full_path = original_folder_loc + filename;
        //readgranule_oneband(original_paths[j+original_lag+smooth_inds[0]],clear_samples,j,"sea_surface_temperature");
        get_var(full_path, original_sst, "sea_surface_temperature");
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                reinstated_clear(y,x,i) = clear_samples(y,x,i+collated_smooth_lag+collated_inds[0]);
                DD = collated_smooth(y,x,i) - original_sst(y,x);
                if(!std::isnan(DD) && DD < T_SMOOTH_COLLATED){
                    reinstated_clear(y,x,i) = original_sst(y,x);
                }
            }
        }
    }
    //SAVENC(reinstated_clear);
    clear_samples.release();
    //approximate using new clear and smooth collated
    //collate one last time
    Mat1f collated_approx(3,dims);
    printf("calculating approximation for brightness_temperature_08um6\n");
    calc_approximate(collated_smooth,reinstated_clear, land_mask,invalid_mask,collated_approx,dims[2]);
    printf("finished approximation for brightness_temperature_08um6\n");
    //Mat1f collated_approx(3,dims);
    //save_mat(collated_approx_paths, collated_approx, "sea_surface_temperature",true);

    collated_paths.clear(); collated_interp_paths.clear();collated_inds.clear();

    printf("starting generation of collated\n");
    //generate save paths for collated_values as well as indicies of hour
    for(i=COLLATED_LAG;i<collated_approx_paths.size() - COLLATED_LAG;i++){
        filename = generate_filename(collated_approx_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            collated_paths.push_back(filename);
            collated_inds.push_back(i);
        }
    }
    
    printf("starting generation of collated interp\n");
    //generate filenames for interp collated
    for(i=collated_inds[0];i<collated_approx_paths.size()-COLLATED_LAG;i++){
        filename = generate_filename(approx_paths[i]); 
        save_loc = collated_folder_loc + filename; 
        printf("approx filename = %s\n",save_loc.c_str());    
        collated_interp_paths.push_back(save_loc);
    }

    collated_size = collated_paths.size();
    collated_interp_size = collated_interp_paths.size();
    dims_collated[2] =collated_size;

    printf("starting collated interp \n");
    interp_inds.clear();
    for(i=0;i<collated_interp_size;i++){
        printf("filename = %s\n",collated_interp_paths[i].c_str());
        filename = generate_filename(collated_interp_paths[i]);
        if(filename[11] == '0' && filename[12] == '0'){
            interp_inds.push_back(i);
        }
    }
    
    printf("starting gamma init\n"); 
    Mat1f collated_final(3,dims_collated);
    collated_final.setTo(NAN);
    Gamma.resize(collated_size,collated_size);
    Gamma.setZero();

    printf("(starting GAMMA\n" );
    //fill gamma matrix
    //first and last row have different values
    Gamma(0,0) = GAMMA;Gamma(0,1) = -GAMMA;
    Gamma(collated_size-1,collated_size-1) = GAMMA; Gamma(collated_size-1,collated_size-2) = -GAMMA;
    for(i=1;i<collated_size-1;i++){
        Gamma(i,i)   = 2*GAMMA;
        Gamma(i,i-1) = -GAMMA;
        Gamma(i,i+1) = -GAMMA;
    }

    printf("Starting Collation on brightness_temperature_08um6\n");
    collate_samples_final(reinstated_clear, collated_approx, collated_final, Gamma, collated_inds,land_mask,invalid_mask,interp,ref_file);
    printf("Finished Collation on brightness_temperature_08um6\n");
    reinstated_clear.release();collated_approx.release();
    dims[2] = collated_interp_size;
    
    SAVENC(collated_final);
    Mat1f collated_interp_final(3,dims);
    collated_interp_final.setTo(NAN);
    for(i=0;i<interp_inds.size();i++){
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                collated_interp_final(y,x,interp_inds[i]) = collated_final(y,x,i);
            }
        }
        
    }
    
    printf("starting interpolation of collated values\n");
    interpolate_hermite(collated_interp_final,land_mask,invalid_mask,collated_interp_size,T_INTERP,2*INTERP_DIST,false);
    remove_last_value(collated_interp_final,invalid_mask,land_mask,collated_interp_size);
    save_mat(collated_interp_paths, collated_interp_final, "sea_surface_temperature",true);
    
    /*
    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "brightness_temperature_10um4",approx, time_size);


    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_approx(save_loc, save_slice,"brightness_temperature_10um4", false);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;

    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "brightness_temperature_11um2",approx, time_size);
    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_approx(save_loc, save_slice,"brightness_temperature_11um2", false);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;

    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "brightness_temperature_12um3",approx, time_size);
    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_approx(save_loc, save_slice,"brightness_temperature_12um3", false);
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;
    approximate(original_paths, smooth_paths,clear_masks, land_mask, invalid_mask, clear_samples, smooth_samples, "sea_surface_temperature",approx, time_size);
    for( i=0;i<time_size-1;i++){
       
        filename = generate_filename(smooth_paths[file_count]);
        printf("filename = %s\n",filename);
        file_count++;

        save_loc = folder_loc+filename;

        //fill_nans_3d(bt11_approx,test_slice,i);
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                save_slice(y,x) = approx(y,x,i);    
            }
        }
        //save_test_nc_fullbands(bt08_slice,bt10_slice,bt11_slice,bt12_slice,sst_slice,save_loc.c_str());
        save_approx(save_loc, save_slice,"sea_surface_temperature", false);
        approx_paths.push_back(save_loc.c_str());
        printf("generated file %s\n", save_loc.c_str());

    }
    file_count = 0;
    */
}

