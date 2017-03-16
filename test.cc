void
test_3dmats()
{
	int y,x,z, i=0;
	int dims[3] = {3,3,3};
	Mat1f test(3,dims);

	for(y=0;y<dims[0];y++){
		for(x=0;x<dims[1];x++){
			for(z=0;z<dims[2];z++){
				test(y,x,z) = i;
				i++;
				printf("value = %d\n",i);
			}
		}
	}

	floatr* p = (float *)test.data;
	size_t elem_step = test.step / sizeof(float);
	for(i = 0; i< 27 ; ++i){
		printf("value = %f\n",++);
	}
}