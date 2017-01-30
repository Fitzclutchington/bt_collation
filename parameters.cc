const int WIDTH = 5500;
const int HEIGHT = 5500;
const int TIMEINT = 9;
const int MAXDIMS = 5;
const bool DEBUG = true;
const float OFFSET = 273.15;
const float SCALE = 0.1;
const short FILLVALUE = 1<<15;
const float PR_CLEAR = 0.1;
const int   MIN_POINTS = 6;
const float DELTA_SST = 0.1;
const float T_SST_SMOOTH = 0.3;

const int FILTER_TIME_SIZE = 3;
const int FILTER_WINDOW_LAG = FILTER_TIME_SIZE/2;
const int SECOND_PASS_SIZE = 25;
const int SECOND_PASS_LAG = SECOND_PASS_SIZE / 2;
const int SMOOTH_TIME_WINDOW = 19;
const int SMOOTH_WINDOW_LAG = SMOOTH_TIME_WINDOW/2;

const float T_GRAD_TIME = 0.05;
const float T_DIAG = 2*DELTA_SST*sqrt(FILTER_TIME_SIZE);
const float T_NN = 2*DELTA_SST;
const float T_BTZ = 0.5;
const float MIN_TEMP = 270;
const float BT08_BT11_DIFF = 4;
const float BT11_BT12_DIFF = 1;
const float T_DT = -5;
const float T_GRAD = 0.3;
const float PASS2 = DELTA_SST;
const float T_INTERP = 0.5;
const float T_COLLATED = 0.5;
const float T_APPROX = 0.3;

const float MU_CLEAR = 0.6;
const float MU_APPROX = 0.4;

const float MU_CLEAR_ABOVE = 0.9;
const float MU_APPROX_ABOVE = 0.1;
const float MU_SMOOTH = 0.1;
const int SUBSIZE = 7;
const int LAND_KERNEL = 7;

const float T_DERIV = DELTA_SST/2;
const float T_EIGEN = 0.1;
const float RATIO = 0.8;
const float T_RATIO = 0;
const float T_COLD = 0.2;
const int SST_RANGE = 3;

const float DT = 1;
const int DL = 4;