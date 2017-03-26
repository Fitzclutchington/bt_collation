//sst_ref, sza, lon- comes from given acspo file 
// lut- (K,L,time) matrix that is the look up table
// K= number of temp bins
// L = number of longitude bins
// ind - what time step

void
generateLUT(const Mat1f &sst_ref,const Mat1f &lon,const Mat1f &sst,Mat1f &lut,
	        const vector<int> ll, const vector<float> temps, const vector<pair<int,int>> valid_points, int ind)
{
	int y,x,k,l,i;
	int K = temps.size();
	int L = ll.size();
	float count, sum;

	for(k=0;k<K;k++){
		for(l=0;l<L;l++){
			count=sum=0;
			for(i=0;i<valid_points.size();i++){
				y = valid_points[i].first;
				x = valid_points[i].second;
				if(lon(y,x)>= (ll[l]-DL) && lon(y,x)<=(ll[l]+DL+1) && sst_ref(y,x)>=(temps[k]-DT) && sst_ref(y,x)<=(temps[k]+DT+1)){
					sum+=clear(y,x);
					count+=1;
				}
			}
			lut(k,l,ind) = NAN;
			if(count>0){
				lut(k,l,ind) = sum/count;
			}
		}
	}
}

// open reference
// get sst_reynolds, sza, lon
// create matrix for lut
// for each time step
//     generateLUT at ind of time_step
// save lut
void
setupLUT(vector<string> &second_pass_paths, vector<string> &original_files, string reference_file){
	
	int time_size = second_pass_paths.size();
	float temp = 270;
	int lon = 59;
	int sza_thresh = 78;
	int y,x,i,j;

	Mat1f sst_ref(HEIGHT,WIDTH);
	Mat1f sza(HEIGHT,WIDTH);
	Mat1f lons(HEIGHT,WIDTH);
	Mat1f sst(HEIGHT,WIDTH);
	Mat1b clear_mask(HEIGHT,WIDTH);
	//Mat1f clear_samples(HEIGHT,WIDTH);

	vector<float> temps;
	vector<int> ll;

	vector<pair<int,int>> valid_points;
	vector<pair<int,int>> sza_points;

	get_var(reference_file.c_str(),sst_ref,"sst_reynolds");
	get_var(reference_file.c_str(),sza,"solar_zenith_angle");
	get_var(reference_file.c_str(),lons,"longitude");

	for(y=0;y<HEIGHT;y++){
		for(x=0;x<WIDTH;x++){
			//only positive longitude
			if(lons(y,x) < 0){
				lons(y,x) += 360;
			}
			//save indicies where sza passed threshold
			if(fabs(sza(y,x) < sza_thresh)){
				sza_points.push_back((y,x));
			}
		}
	}

	while(temp <= 310){
		temps.push_back(temp);
		temp+=DT;
	}

	while(lon<=222){
		ll.push_back(lon);
		lon+=DL;
	}

	int dims[3] = {temps.size(),ll.size(),time_size};
	Mat1f lut(3,dims);
	lut.setTo(NAN);
	printf("computed values for ll and temps\n");
	for(i=0;i<time_size;i++){
		read_mask(second_pass_paths[i],clear_mask,-1);
    	get_var(original_files[i+FILTER_WINDOW_LAG+SECOND_PASS_LAG],sst, "sea_surface_temperature");
    	//apply_mask_2d(clear_mask, sst);
    	valid_points.clear();
    	for(j=0;j<sza_points.size();j++){
    		if(!std::isnan(clear_mask(sza_points[j].first,sza_points[j].second))){
    			valid_points.push_back(sza_points[j]);
    		}
    	}
    	generateLUT(sst_ref,lons,sst,lut,ll, temps, valid_points, i);
    	printf("finished time %d\n",i);
	}
	SAVENC(lut);
}