
void
smoothing3D_mat(const Mat1f &sst, const Mat1f &reference, const Mat1b &land_mask, const Mat1b &invalid_mask, Mat1f &sst_3dsmooth, int *window)
{
	int x,y,t,i,j,k;

	int threshold;

	threshold = PR_CLEAR * (2*window[0]+1)*(2*window[1]+1)*(2*window[2]+1);
	//threshold =1;
	//printf("%d \n",threshold);
	printf("starting loop\n");

	//windowed_nanmean_3d(bt08, land_mask, invalid_mask, bt08_3dsmooth, window, threshold);
	//windowed_nanmean_3d(bt10, land_mask, invalid_mask, bt10_3dsmooth, window, threshold);
	//windowed_nanmean_3d(bt11, land_mask, invalid_mask, bt11_3dsmooth, window, threshold);
	//windowed_nanmean_3d(bt12, land_mask, invalid_mask, bt12_3dsmooth, window, threshold);
    //subtract ref from sst ( sst - ref)
    
	windowed_nanmean_3d(sst, land_mask, invalid_mask,  sst_3dsmooth, window, threshold);
    //add back reference to sst_3dsmooth ( sst_3dsmooth + ref)
    for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            sst_3dsmooth(y,x) = sst_3dsmooth(y,x) + reference(y,x);
        }
    }
    
}



void
smooth_samples_collated(const Mat1f &samples, Mat1f &smooth_samples,const Mat1b &land_mask, const Mat1b &invalid_mask,
                        const string ref_file, int *window, int time){
	
	int j,y,x,t;
	int i=0;
    int time_size = window[2]*2 + 1;
    printf("time size = %d\n",time_size);
    int time_ind = time_size;
    int file_count = 0;
    int start_ind = 0;
    int dims[3] = {HEIGHT,WIDTH,time_size};
    int sample_size = time;
    string folder_loc = "data/smooth";
 	string filename, save_loc;


    Mat1f sst(3,dims); 

    Mat1b clear_mask(HEIGHT,WIDTH);


  	Mat1f sst_3dsmooth(HEIGHT, WIDTH);
    Mat1f reference(HEIGHT,WIDTH);


    read_reference(ref_file,reference,"sst_reynolds");
 
	for(y=0;y<HEIGHT;y++){
        for(x=0;x<WIDTH;x++){
            for(t=0;t<time_size;t++){
                sst(y,x,t) = samples(y,x,t) - reference(y,x);
            }
        }
    }
    


    //SAVENC(bt08);
    printf("read granules\n");
    while(time_ind < sample_size+1){     
 
        printf("starting 3d smoothing \n");

        sst_3dsmooth.setTo(NAN);
        smoothing3D_mat( sst, reference, land_mask, invalid_mask, sst_3dsmooth, window);
        printf("finished 3d smoothing\n");
        
        for(y=0;y<HEIGHT;y++){
            for(x=0;x<WIDTH;x++){
                smooth_samples(y,x,file_count) = sst_3dsmooth(y,x);
            }
        }  

    	file_count++;

        
        if(time_ind<sample_size){

            for(y=0;y<HEIGHT;y++){
                for(x=0;x<WIDTH;x++){
                    sst(y,x,start_ind) = samples(y,x,time_ind) - reference(y,x);                    
                }
            }

	        printf("read granule %d\n",time_ind);        
	        start_ind = (start_ind +1) % time_size;        

	      }
	    time_ind++; i++;
    }
        
}