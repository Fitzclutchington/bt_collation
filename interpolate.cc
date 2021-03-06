void
interpolate(Mat1f &smooth,  const Mat1b &land_mask, const Mat1b &invalid_mask, int time_size, float threshold_y, float threshold_x, bool thresh, bool window)
{
	int x,y,i,t,j;
	int w = 1;
	int xi_size = 0; int xd_size = 0;
	float val=0,count = 0;
	vector<int> xd;
	vector<float> yd;
	vector<int> xi;
	vector<float> yi;

	for(y=0;y<HEIGHT;y++){
		for(x=0;x<WIDTH;x++){
			if(invalid_mask(y,x) == 0 && land_mask(y,x) == 0){
				xd.clear();
				yd.clear();
				xi.clear();
				yi.clear();

				//Fill xd,xi,yd
				for(t=0;t<time_size;t++){
					
					if(std::isnan(smooth(y,x,t))){
						xi.push_back(t);
					}
					else{
						if(thresh){
							if((!std::isnan(smooth(y,x,t+1)) && fabs(smooth(y,x,t) - smooth(y,x,t+1)) < T_DERIV) || std::isnan(smooth(y,x,t+1))){
								xd.push_back(t);
								if(window){
									for(i=-w;i<w+1;i++){
										for(j=-w;j<w+1;j++){
											if(!std::isnan(smooth(y+i,x+j,t))){
												val+=smooth(y+i,x+j,t); 
												count++;
											}
										}
										yd.push_back(val/count);
										val =0;
										count =0;
									}
								}
								else{
									yd.push_back(smooth(y,x,t));
								}
							}
							else{
								smooth(y,x,t) = NAN;
								xi.push_back(t);
							}
						}
						else{
							xd.push_back(t);
							if(window){
								for(i=-w;i<w+1;i++){
									for(j=-w;j<w+1;j++){
										if(!std::isnan(smooth(y+i,x+j,t))){
											val+=smooth(y+i,x+j,t); 
											count++;
										}
									}
									yd.push_back(val/count);
									val =0;
									count =0;
								}
							}
							else{
								yd.push_back(smooth(y,x,t));
							}
						}
					}
					
				}
				
				xi_size = xi.size();
				xd_size = xd.size();

				if(xd_size && xi_size){
					pwl_interp_1d(xd,yd,xi,yi,threshold_y,threshold_x);
					for(i=0;i<xi_size;i++){
						smooth(y,x,xi[i]) = yi[i];
					}
				}

		
			}
		}
	}
}

void
interpolate_hermite(Mat1f &smooth,  const Mat1b &l2p_mask, int time_size, float y_threshold, float x_threshold, bool thresh)
{
	int x,y,i,t;
	int xi_size = 0; int xd_size = 0;
	vector<double> xd;  // data points
	vector<double> yd; //data values
	vector<double> xi; //interpolation points
	vector<double> yi; // interpolation values

	for(y=0;y<HEIGHT;y++){
		for(x=0;x<WIDTH;x++){
			if(l2p_mask(y,x) == 0){
				xd.clear();
				yd.clear();
				xi.clear();
				yi.clear();

				//Fill xd,xi,yd
				for(t=0;t<time_size;t++){
					
					// first check if smooth value is nan
					if(std::isnan(smooth(y,x,t))){
						xi.push_back(t);
					}

					else{
						// if threshold is true p1 - p2 must be less than T_derive 
						if(thresh){
							if((!std::isnan(smooth(y,x,t+1)) && fabs(smooth(y,x,t) - smooth(y,x,t+1)) < T_DERIV) ){
								xd.push_back(t);
								yd.push_back(smooth(y,x,t));
							}
							else{
								smooth(y,x,t) = NAN;
								xi.push_back(t);
							}
						}
						else{
							xd.push_back(t);
						    yd.push_back(smooth(y,x,t));							
						}
					}
					
				}
				
				xd_size = xd.size();
				xi_size = xi.size();

				if(xd_size > 0 && xi_size >0){
					MonotCubicInterpolator MCI(xd ,yd);
					
					for(i=0;i<xi_size;i++){

						smooth(y,x,xi[i]) = MCI.evaluate(xi[i], y_threshold, x_threshold);
					}
				}

		
			}
		}
	}
}

