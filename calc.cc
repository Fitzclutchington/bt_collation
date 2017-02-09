// Compute the gradient magnitude of image src into dst.
//
// src -- source image
// dst -- destination image of gradient magnitude (output)
// dX, dY -- derivative in the X and Y directions (output)
//
void
gradientmag(const Mat1f &src, Mat1f &dst)
{
    //printf("about to seg fault\n");
    Mat1f dX(HEIGHT,WIDTH);
    Mat1f dY(HEIGHT,WIDTH);
    int y,x;
    Mat h = (Mat_<float>(5,1) <<
            0.036420, 0.248972, 0.429217, 0.248972, 0.036420);
    Mat hp = (Mat_<float>(5,1) <<
              0.108415, 0.280353, 0, -0.280353, -0.108415);

    sepFilter2D(src, dX, -1, h, hp);
    // We negate h here to fix the sign of Dy
    sepFilter2D(src, dY, -1, hp, -h);
    for(y=0;y<HEIGHT;y++){
      for(x=0;x<WIDTH;x++){
          dst(y,x) = sqrt(dX(y,x)*dX(y,x) + dY(y,x)*dY(y,x));
      }
    }
}

void
compute_gradient(const Mat1f &bt11, Mat1f &mags, int *dims, int *inds)
{
    int a,b,i,x,y;
	  Mat1f tempmags(HEIGHT,WIDTH);
    Mat1f tempmat(HEIGHT,WIDTH);
    #pragma omp for
    for( i = 0; i < dims[2];i++){
    	for(a=0;a<HEIGHT;a++){
            for(b=0;b<WIDTH;b++){
                tempmat(a,b) = bt11(a,b,inds[i]);
            }
        }

    	gradientmag(tempmat,tempmags);
        for( y = 0; y < HEIGHT; y++){
        	for( x =0; x < WIDTH; x++){
        		mags(y,x,inds[i]) = tempmags(y,x);
        	}
        }
    }
}

float
compute_nanmean(const Mat1f &data, const int *loc, const int *windows, int *inds)
{
  float sum = 0;
  float d;
    int i,j,k;
  int count = 0;
  int y = loc[0]; int x = loc[1]; int z = loc[2];
  int w1 = windows[0]; int w2 = windows[1]; int w3 = windows[2];

    for(i=-w1;i<w1+1;i++){
        for(j=-w2;j<w2+1;j++){
            for(k=-w3;k<w3+1;k++){
                d = data(y+i,x+j,inds[z+k]);
                if(!std::isnan(d)){
                    count += 1;
                    sum += d;
                }
            }
        }
    }
    return sum/count;
}

void
fill_nans(const Mat1f &bt11_masked, Mat1s &test_slice){
    #pragma omp for
    for(int y = 0;y<HEIGHT;y++){
		for(int x = 0;x<WIDTH;x++){
        	if(std::isnan(bt11_masked(y,x))){
        		test_slice(y,x) = FILLVALUE;
        	}
            else{
            	test_slice(y,x) = round((bt11_masked(y,x)-OFFSET)/SCALE);
        	}
        }
    }
}

void
fill_nans_3d(const Mat1f &bt11_3d, Mat1s &test_slice, int slice){
    #pragma omp for
    for(int y = 0;y<HEIGHT;y++){
    for(int x = 0;x<WIDTH;x++){
          if(std::isnan(bt11_3d(y,x,slice))){
            test_slice(y,x) = (short) FILLVALUE;
          }
            else{
              test_slice(y,x) = (short) (round((bt11_3d(y,x,slice)-OFFSET)/SCALE));
          }
        }
    }
}

void
calculate_diffs(const Mat1f &a, const Mat1f &b, Mat1f &diffs,int cur_ind, bool absval){
  //printf("calculating diffs\n");

  //printf("suxxess at time %d\n",t);
  for(int y = 0;y<HEIGHT;y++){
    for(int x = 0;x<WIDTH;x++){
      diffs(y,x) = NAN;
      if(!std::isnan(a(y,x,cur_ind)) && !std::isnan(b(y,x,cur_ind))){
        if(absval){
            diffs(y,x) = fabs(a(y,x,cur_ind) - b(y,x,cur_ind));
        }
        else{
            diffs(y,x) = a(y,x,cur_ind) - b(y,x,cur_ind);
        }
      }
    }
  }
}   


void
calculate_sums(const Mat1f &a, const Mat1f &b, Mat1f &sums,int* dims, int *inds){
  //printf("calculating diffs\n");
  int time = dims[2];
    for(int t = 0; t<time; t++){
      //printf("suxxess at time %d\n",t);
      for(int y = 0;y<HEIGHT;y++){
        for(int x = 0;x<WIDTH;x++){
                sums(y,x,inds[t]) = a(y,x,inds[t]) + b(y,x,inds[t]);            
          }
        }
    }
    
}

void 
pwl_interp_1d (const vector<int> &xd, const vector<float> &yd, const vector<int> &xi, 
                vector<float> &yi, float threshold)


//  Purpose:
//
//    PWL_INTERP_1D evaluates the piecewise linear interpolant.
//
//  Discussion:
//
//    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
//    linear function which interpolates the data (XD(I),YD(I)) for I = 1
//    to ND.
//
//  Parameters:
//
//    Input, vector<int> xd, the data points.
//
//    Input, vector<float> yd, the data values.
//
//    Input, vector<float> yr_r, the data points of reference
//
//    Input, vector<int> xi, the interpolation points.
//
//
{
  int i;
  int k;
  int pos =0;
  float t;
  int nd = xd.size();
  int ni = xi.size();
  //printf("nd = %d ni = %d\n",nd,ni);
  //if there is only one data point, set all interpolation
  //points to the value at that point
  if ( nd == 1 ){
    for ( i = 0; i < ni; i++ ){
      yi.push_back(NAN);
    }
  }

  else{
    //for xi before xd interval
    while(xi[pos]<xd[0] && pos < ni){
        yi.push_back(NAN);
        pos++;
    }

    if(xi[pos] < xd[nd-1] && pos < ni){
      for ( k = 0; k< nd-1; k++ ){

        if(xi[pos] < xd[k+1]){
          if(fabs(yd[k+1] - yd[k]) < threshold){

            while(xi[pos] < xd[k+1]){

              t = ( xi[pos] - xd[k+1] ) / (float)( xd[k] - xd[k+1] );
              yi.push_back(( 1.0 - t ) * yd[k+1] + t * yd[k]);
              pos++;
            }
          }
          else{

            while(xi[pos] < xd[k+1]){
              
              yi.push_back(NAN);
              pos++;
            }
          }
        }
      }
    }

    while(xi[pos] > xd[nd-1] && pos < ni){
        //printf("outside loop pos=%d\n",pos);
        yi.push_back(NAN);
        pos++; 
    }
  }
}


void
compute_indicies(int *indicies, int start_ind, int time_size){
  int i;
  int val;
  for(i=0;i<time_size;i++){
    val = ((i-start_ind)%time_size + time_size) % time_size;
    indicies[i] = val;
  }
}

string
convert_int_to_string(int val){
  std::stringstream ss; 
  ss << val;
  string str = ss.str();
  return str;
}

void
windowed_nanmean_3d(const Mat1f &samples, const Mat1b &land_mask, const Mat1b &invalid_mask, const Mat1b &ice_masks, 
                    Mat1f &nanmean, int *window, int *inds, float threshold)
{
  Mat1f time_sum(HEIGHT,WIDTH);
  Mat1f time_count(HEIGHT,WIDTH);

  int y,x,t,i,j;
  int y_dim = window[0];
  int x_dim = window[1];
  int t_dim = window[2];

  int t_len = 2*t_dim +1;

  float count,sum,left_sum,left_count;

  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      for(t=0;t<t_len;t++){
        if(!std::isnan(samples(y,x,inds[t]))){
          time_sum(y,x) += samples(y,x,inds[t]);
          time_count(y,x) += 1;
        }
      }
    }
  }
  printf("finished nan sum of time\n");

  //calculate stats in time dimension
  for(y=y_dim;y<HEIGHT-y_dim;y++){
    x=x_dim;
    sum = 0;
    count  = 0;
    left_sum = 0;
    left_count = 0;
    //find first first valid pixel
    while((land_mask(y,x) !=0 || invalid_mask(y,x) != 0) ){
    // && ice_masks(y,x,t_dim) != 0){
      x++;
    }
        // calc first window sum and count

    for(i=-y_dim;i<y_dim+1;i++){
      for(j=-x_dim;j<x_dim+1;j++){
        //if(y+i > HEIGHT || x +j > WIDTH){
        //  printf("x = %d , y= %d\n", y+i, x+j);
        //}          
        count+= time_count(y+i,x+j);
        sum += time_sum(y+i,x+j);         
        
      }
    }


    if(count > threshold && x<WIDTH-x_dim){

      nanmean(y,x) = sum /count;
    }

    // now repeat for all other valid x values
    for(x=x+1;x<WIDTH-x_dim;x++){
      //remove all values of previous left
      sum -= left_sum;
      count -= left_count;

      //calculate new left and new right
      left_sum = 0;
      left_count = 0;
      for(i=-y_dim;i<y_dim+1;i++){  
        left_count+=time_count(y+i,x-x_dim);
        left_sum += time_sum(y+i,x-x_dim);  
        count+= time_count(y+i,x+x_dim);
        sum += time_sum(y+i,x+x_dim);     
      }


      //printf("count = %f sum = %f\n", count, sum);
      if(land_mask(y,x)==0 && invalid_mask(y,x) ==0 && ice_masks(y,x,t_dim) ==0){
        // calculate this pixels value at 3dsmooth(y,x)
        if(count>threshold){
          nanmean(y,x) = sum /count;
        }
      }
    }
  }
}

void
windowed_nanmean_2nd_pass(const Mat1f &samples, const Mat1f &counts, const Mat1b &land_mask, const Mat1b &invalid_mask, 
                          Mat1f &nanmean, int *window, int mid)
{
  Mat1f time_sum(HEIGHT,WIDTH);
  Mat1f time_count(HEIGHT,WIDTH);

  time_sum.setTo(0);
  time_count.setTo(0);
 
  int y,x,t,i,j,time;
  int y_dim = window[0];
  int x_dim = window[1];
  int t_dim = window[2];

  float count,sum,left_sum,left_count;

  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      for(t=-t_dim;t<t_dim+1;t++){
        time =(((mid + t) % SECOND_PASS_SIZE + SECOND_PASS_SIZE)%SECOND_PASS_SIZE);
        if(!std::isnan(samples(y,x,time))){
          time_sum(y,x) += samples(y,x,time);
          time_count(y,x) += 1;
        }
      }
    }
  }

  printf("finished nan sum of time\n");
  //calculate stats in time dimension
  for(y=y_dim;y<HEIGHT-y_dim;y++){
    x=x_dim;
    sum = 0;
    count  = 0;
    left_sum = 0;
    left_count = 0;
    //find first first valid pixel
    while(land_mask(y,x) !=0 && invalid_mask(y,x) != 0){
      x++;
    }
    //printf("first valid x = %d\n",x);
    // calc first window sum and count    
    for(i=-y_dim;i<y_dim+1;i++){
      for(j=-x_dim;j<x_dim+1;j++){          
        count+= time_count(y+i,x+j);
        sum += time_sum(y+i,x+j);         
        
      }
    }

    if(counts(y,x) >0){
      nanmean(y,x) = sum /count;
    }
    

    // now repeat for all other valid x values
    for(x=x+1;x<WIDTH;x++){
      //remove all values of previous left
      sum -= left_sum;
      count -= left_count;

      //calculate new left and new right
      left_sum = 0;
      left_count = 0;
      for(i=-y_dim;i<y_dim+1;i++){
        left_count+=time_count(y+i,x-x_dim);
        left_sum += time_sum(y+i,x-x_dim);  
        count+= time_count(y+i,x+x_dim+1);
        sum += time_sum(y+i,x+x_dim+1);     
      }


      //printf("count = %f sum = %f\n", count, sum);
      if(land_mask(y,x)==0 && invalid_mask(y,x) ==0 && counts(y,x) > 0){
        // calculate this pixels value at 3dsmooth(y,x)
        nanmean(y,x) = sum /count;
      }
    }
  }
}

void
get_landborders(const Mat1b &land_mask, Mat1b &border_mask,int kernel_size)
{
  int x,y;
  Mat element = getStructuringElement( MORPH_RECT,
                                       Size( kernel_size, kernel_size ) );
  dilate(land_mask, border_mask,element);
  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      border_mask(y,x) -= land_mask(y,x);
    }
  }
}


void
compute_eigenvals(const Mat1f &bt08,const Mat1f &bt10,const Mat1f &bt11,const Mat1f &bt12,
                  const Mat1b border_mask, Mat1f &clear_samples,int cur_ind)
{
  int y,x,t,i,j,k;
  int y_delta = 2;
  int x_delta = 2;
  int t_delta = 2;

  int min_num = (2*y_delta +1) *(2*x_delta + 1)*(2*t_delta+1)/2;
  float bt08_sum,bt10_sum,bt11_sum,bt12_sum,count,window_sum,row_sum, res_mean;
  float temp_bt08;
  float temp_bt10;
  float temp_bt11;
  float temp_bt12;

  float bt08_mean;
  float bt10_mean;
  float bt11_mean;
  float bt12_mean;

  Mat1f counts(HEIGHT,WIDTH);
  Mat1f vals(HEIGHT,WIDTH);

  Mat1b eigen_mask(HEIGHT, WIDTH);
  eigen_mask.setTo(0);
  vector<float> valid_bt08;
  vector<float> valid_bt10;
  vector<float> valid_bt11;
  vector<float> valid_bt12;

  vector<int> left_inds;

  Vector4f ones(1,1,1,1);
  Vector4f e1;
  MatrixXf r;
  Matrix4f A;

  for(y=y_delta;y<HEIGHT-y_delta;y++){
    for(x=x_delta;x<WIDTH-x_delta;x++){
      if(!std::isnan(clear_samples(y,x)) && border_mask(y,x) == 0){
    
        // calc first window
        // we know that first left are nans so we don't calculate left inds     
        bt08_sum=bt10_sum=bt11_sum=bt12_sum=0;
        valid_bt08.clear();
        valid_bt10.clear();
        valid_bt11.clear();
        valid_bt12.clear();
        for(i=-y_delta;i<y_delta+1;i++){
          for(j=-x_delta;j<x_delta+1;j++){              
            for(k=-t_delta;k<t_delta+1;k++){
              t = ((((cur_ind+k)%FILTER_TIME_SIZE)+FILTER_TIME_SIZE) % FILTER_TIME_SIZE);
              temp_bt08 = bt08(y+i,x+j,t);
              temp_bt10 = bt10(y+i,x+j,t);
              temp_bt11 = bt11(y+i,x+j,t);
              temp_bt12 = bt12(y+i,x+j,t);
              if(!std::isnan(temp_bt08) && !std::isnan(temp_bt10) && !std::isnan(temp_bt11) && !std::isnan(temp_bt12)){
                valid_bt08.push_back(temp_bt08);
                valid_bt10.push_back(temp_bt10);
                valid_bt11.push_back(temp_bt11);
                valid_bt12.push_back(temp_bt12);

                bt08_sum+= temp_bt08;
                bt10_sum+=temp_bt10;
                bt11_sum+=temp_bt11;
                bt12_sum+=temp_bt12;
              }
            }
          }
        }
  
        //if numberof pixels in window is greater tan threshold
        // calculate the mean of the norm of the pixels
        // projected into the second eigenvector
        count = valid_bt08.size();
        counts(y,x) = count;
        if(valid_bt08.size() > min_num){
          bt08_mean =bt08_sum/count;
          bt10_mean =bt10_sum/count;
          bt11_mean =bt11_sum/count;
          bt12_mean =bt12_sum/count;

          MatrixXf window(valid_bt08.size(),4);
          for(i=0;i<valid_bt08.size();i++){
            window(i,0) = valid_bt08[i] - bt08_mean;
            window(i,1) = valid_bt10[i] - bt10_mean;
            window(i,2) = valid_bt11[i] - bt11_mean;
            window(i,3) = valid_bt12[i] - bt12_mean;
          }
          
          A = (window.transpose()*window);
          e1 = A*(A*ones);
          e1/=sqrt(e1.transpose()*e1);
          r = window - (window*e1)*e1.transpose();
          window_sum =0;
          for(i=0;i<valid_bt08.size();i++){
            row_sum = 0;
            row_sum+=r(i,0)*r(i,0);
            row_sum+=r(i,1)*r(i,1);
            row_sum+=r(i,2)*r(i,2);
            row_sum+=r(i,3)*r(i,3);
            row_sum = sqrt(row_sum);
            window_sum += row_sum;
          }
          
          res_mean = window_sum/(float)valid_bt08.size();
          vals(y,x) = res_mean;
          if(res_mean > T_EIGEN){
            clear_samples(y,x) = NAN;
          }
          else{
            eigen_mask(y,x) = 1;
          }
        }
      }
    }
  }
  //SAVENC(eigen_mask);
}


void
calculate_bt_ratio(const Mat1f &bt08, const Mat1f &bt11,const Mat1f &bt12, 
                   Mat1f &bt_ratio,float ratio,int ind)
{
  int y,x;
  for(y=0;y<HEIGHT;y++){
    for(x=0;x<WIDTH;x++){
      bt_ratio(y,x) = bt08(y,x,ind) + 0.8*bt11(y,x,ind) - (1+0.8)*bt12(y,x,ind);
    }
  }
}

string
generate_filename(const string file_loc)
{
  string div = "/";
  int temp_pos = 0;
  int pos = temp_pos;
  string filename;
  while(temp_pos != string::npos){
    temp_pos = file_loc.find(div,pos+1);
    if(temp_pos != string::npos){
      pos = temp_pos;
    }
  }

  filename = file_loc.substr(pos);
  return filename;
}
