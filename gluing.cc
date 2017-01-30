/*
void
combine_clear_approx(const vector<string> clear_paths, const vector<string> approx_paths)
{
	Mat1f glued(HEIGHT,WIDTH);
	Mat1f clear(HEIGHT,WIDTH);
	Mat1f approx(HEIGHT,WIDTH);

	Mat1b flag(HEIGHT,WIDTH);
	string filename;

	for(i=0;i<smooth_paths.size();i++){
		readgranule_oneband(approx_paths[i],smooth);
		readgranule_oneband(clear_paths[i+10],clear);

		for(y=0;y<HEIGHT;y++){
			for(x=0;x<WIDTH;x++){
				if(!(std::isnan(clear(y,x)))){
					glued(y,x)=clear(y,x);
					flag(y,x) = 0;
				}
				else{
					glued(y,x)=approx(y,x);
					flag(y,x) = 255;
				}
			}
		}

		filename = "data/glued/glued_samples" + convert_int_to_string(i) +".nc";
        
        save_test_nc_final(glued,flag,filename.c_str());
        printf("generated file %s\n", filename.c_str());
	}
}
*/