void
apply_l2p_flags(const Mat1b &ice_mask,const Mat1b &land_mask,const Mat1b &invalid_mask,Mat1f &bt11_clear, int ind)
{
    int x,y;
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(ice_mask(y,x,ind)==255 || land_mask(y,x)==255 || invalid_mask(y,x)==255){
                bt11_clear(y,x) = NAN;
            }
        }
    }
}

/*

Cloud mask determined by the difference of neighboring pixels diagonal

*/
void 
compute_diagmask(const Mat1f &samples, Mat1f &clear_samples, int cur_ind)
{
    int x,y,i, ind;
    int w = FILTER_WINDOW_LAG;
    float DD,prod;
    Mat1b diag_mask(HEIGHT,WIDTH);
    diag_mask.setTo(0);
    //Mat1f DD_test(3,dims);
    Mat1f window(1,FILTER_TIME_SIZE);
    Mat1f u = Mat1f::ones(1,FILTER_TIME_SIZE);
    Mat1f powers;
    Mat1f D(1,FILTER_TIME_SIZE);
    u = u/sqrt(sum(u).val[0]);
    //printf("starting diag loop\n");
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(!std::isnan(clear_samples(y,x))){
    
                prod = 0;
                DD=0;
                for(i=0;i<FILTER_TIME_SIZE;i++){
                    ind = ((cur_ind-w+i % FILTER_TIME_SIZE) + FILTER_TIME_SIZE) % FILTER_TIME_SIZE;
                    window(i) = samples(y,x,ind);
                    prod+=(window(i)*u(i));
                }
                
                for(i=0;i<FILTER_TIME_SIZE;i++){
                    D(i) = window(i) - prod*u(i);
                    D(i) = D(i)*D(i);
                    DD+=D(i);
                }
                if(DD > T_DIAG || std::isnan(DD)){
                    clear_samples(y,x) = NAN;
                }
                else{
                    diag_mask(y,x) = 1;
                }
            }
        }
    }
    //SAVENC(diag_mask);
}


/*

Cloud mask determined by the difference of neighboring pixels

*/
void 
compute_nnmask(const Mat1f &samples,Mat1f &clear_samples, int cur_ind)
{
    int x,y, ind_next, ind_prev;
    float D_right, D_left;  
    ind_next = (cur_ind+1) % FILTER_TIME_SIZE; 
    ind_prev = ((cur_ind-1 % FILTER_TIME_SIZE) + FILTER_TIME_SIZE) % FILTER_TIME_SIZE;  
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(!std::isnan(clear_samples(y,x))){                
                D_right = fabs(samples(y,x,cur_ind) - samples(y,x,ind_next));
                D_left = fabs(samples(y,x,cur_ind) - samples(y,x,ind_prev));
                if(D_right > T_NN || D_left > T_NN || (std::isnan(D_right) && std::isnan(D_left))){
                    clear_samples(y,x) = NAN;
                }
            }
        }
    }
    //SAVENC(nn_mask);
}

/*

Cloud mask determined by a pixels temperature
255 = ocean
*/
void 
compute_threshmask(const Mat1f &samples, Mat1f &samples_clear, int cur_ind,  float threshold,bool direction)
{
    int x,y;
    float DD;    
    //printf("success at time %d",t);
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){   

            DD = samples(y,x,cur_ind);
            
                if(direction){
                    if(DD <threshold){
                        samples_clear(y,x) = NAN;
                    }
                }
                else{
                    if(DD >threshold){
                        samples_clear(y,x) = NAN;
                    }
                } 
            }               
        
    }
}


//posibly make diffs mask instead
void 
compute_threshmask_2d(const Mat1f &samples, Mat1f &samples_clear, int cur_ind,  float threshold,bool direction)
{
    int x,y;
    float DD;    
    //printf("success at time %d",t);
    Mat1b ratio_mask(HEIGHT,WIDTH);
    ratio_mask.setTo(0);
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){               
            DD = samples(y,x);
            if(!std::isnan(DD)){
                if(direction){
                    if(DD < threshold){
                        samples_clear(y,x) = NAN;
                    }
                    else{
                        ratio_mask(y,x) = 1;
                    }
                }
                else{
                    if(DD > threshold){
                        samples_clear(y,x) = NAN;
                    }
                    else{
                        ratio_mask(y,x) = 1;
                    }
                }                
            }
        }
    }
    //SAVENC(ratio_mask);
}


void
compute_btz_mask(const Mat1f &bt08, const Mat1f &bt10, const Mat1f &bt11, const Mat1f &bt12, float *btweights, Mat1f &bt11_clear, int cur_ind)
{
  int x,y;
  float btz;
  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      btz = btweights[0]*bt08(y,x,cur_ind) + btweights[1]*bt10(y,x,cur_ind)+ btweights[2]*bt11(y,x,cur_ind) + btweights[3]*bt12(y,x,cur_ind);
      if(btz>T_BTZ){
        bt11_clear(y,x) = NAN;
      }
    }
  }
}

void
compute_cold_mask(const Mat1f &bt08, const Mat1f &bt10, Mat1f &bt11_clear, int cur_ind)
{
  int x,y;
  float cold;
  Mat1b cold_mask(HEIGHT,WIDTH);
  cold_mask.setTo(0);
  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      cold = (100*(bt10(y,x,cur_ind) - bt08(y,x,cur_ind)))/((bt10(y,x,cur_ind) + bt08(y,x,cur_ind)));
      if(cold<T_COLD){
        bt11_clear(y,x) = NAN;
      }
      else{
        cold_mask(y,x) = 1;
      }
    }
  }
  //SAVENC(cold_mask);
}


void
filter_histogram(const Mat1f &bt08, const Mat1f &bt10, const Mat1f &bt11, const Mat1f &bt12, const Mat1f &lats, Mat1f &bt11_clear, int cur_ind)
{
    int x,y;
    int r1_L = -3;   int r1_R = 7;
    int r3_L = 1;    int r3_R = 6;
    int t_L  = 265;  int t_R  = 301;

    float t, l, r1, r3;
    float delta_r = 0.1;
    float delta = 0.5;
    float delta_l = 5;

    float dt_L = t_L/delta; float dt_R = t_R/delta;
    float dr1_L = r1_L/delta_r; float dr1_R = r1_R/delta_r;
    float dr3_L = r3_L/delta_r; float dr3_R = r3_R/delta_r;
    float dlat_L =0; 

    //int dims[4] = {67, 101, 81, 51};
    int dims[4] = {73 , 19 , 101 , 51};
    int d[4] = {0,0,0,0};
    string filename = "LUT.nc";

    Mat1b H4D(4,dims);
    Mat1b test(HEIGHT,WIDTH);
    test.setTo(0);

    open_LUT(filename, H4D, dims);
    
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            t = round(bt11(y,x,cur_ind)/delta);
            if(t>=dt_L && t<=dt_R){
                //printf("passed first if statement\n");
                l = round(fabs(lats(y,x))/delta_l);
                r1=round((bt08(y,x,cur_ind)-bt12(y,x,cur_ind))/delta_r);
                r3=round((bt11(y,x,cur_ind)-bt12(y,x,cur_ind))/delta_r);

                
                if(r1>=dr1_L && r1<=dr1_R && r3>=dr3_L && r3<=dr3_R){
                    //printf("passed second if statement\n");
                    // for latitude
                    d[0] = t-dt_L;
                    d[1] = l-dlat_L;
                    d[2] = r1-dr1_L;
                    d[3] = r3-dr3_L;
                    /* withou lat
                    d[0] = t-dt_L;
                    d[1] = r1-dr1_L;
                    d[2] = r2-dr2_L;
                    d[3] = r3-dr3_L;
                    */
                    test(y,x) = H4D(d);
                    if(!(H4D(d))){
                        //printf("passed third if statement\n");
                        
                        bt11_clear(y,x) = NAN;
                    }
                }
                else{
                    bt11_clear(y,x) = NAN;
                }
            }
            else{
                bt11_clear(y,x) = NAN;
            }
        }
    }

}

void
generate_mask(const Mat1f &clear_samples, Mat1b &clear_mask)
{
    int y,x;
    clear_mask.setTo(0);
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(!std::isnan(clear_samples(y,x))){
                clear_mask(y,x) = 255;
            }
        }
    }
}

void
apply_mask(const Mat1b &mask, Mat1f &bt08, Mat1f &bt10, Mat1f &bt11,Mat1f &bt12, Mat1f &sst, int cur_ind)
{
    int y,x;
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(mask(y,x) == 0){
                bt08(y,x,cur_ind) = NAN;
                bt10(y,x,cur_ind) = NAN;
                bt11(y,x,cur_ind) = NAN;
                bt12(y,x,cur_ind) = NAN;
                sst(y,x, cur_ind) = NAN;
            }
        }
    }
}

void
apply_mask_slice(const Mat1b &mask, Mat1f &samples, int cur_ind, bool dim)
{
    // if dim is true mask is a 3d array
    // if dim is false mask is a 2d array
    // necessary to deal with 3rd index
    int y,x;
    if(dim){
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                if(mask(y,x,cur_ind) == 0){
                    samples(y,x,cur_ind) = NAN;
                }
            }
        }
    }
    else{
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                if(mask(y,x) == 0){
                    samples(y,x,cur_ind) = NAN;
                }
            }
        }
    }
}

void
apply_mask_2d(const Mat1b &mask, Mat1f &samples)
{
    int y,x;
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            if(mask(y,x) == 0){
                samples(y,x) = NAN;
            }
        }
    }
}
