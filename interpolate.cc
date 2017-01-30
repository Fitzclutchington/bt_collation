void
interpolate(Mat1f &smooth,  const Mat1b &land_mask, const Mat1b &invalid_mask, int time_size, float threshold, bool thresh, bool window)
{
	int x,y,i,t,j;
	int w = 1;
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
			
				if(xd.size() > 0 && xi.size()>0){
					pwl_interp_1d(xd,yd,xi,yi,threshold);
					for(i=0;i<xi.size();i++){
						smooth(y,x,xi[i]) = yi[i];
					}
				}

		
			}
		}
	}
}