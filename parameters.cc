const int WIDTH = 5500;
const int HEIGHT = 5500;
const int TIMEINT = 144;
const int MAXDIMS = 5;

const bool DEBUG = true;
const float PR_CLEAR = 0.1;
const int   MIN_POINTS = 6;
const float DELTA_SST = 0.1;
const float T_SST_SMOOTH = 0.3;

const int DIAG_SIZE = 3;
const int DIAG_LAG = DIAG_SIZE/2;

const int FILTER_TIME_SIZE = 5;
const int FILTER_WINDOW_LAG = FILTER_TIME_SIZE/2;
const int SECOND_PASS_SIZE = 25;
const int SECOND_PASS_LAG = SECOND_PASS_SIZE / 2;
const int SMOOTH_TIME_WINDOW = 19;
const int SMOOTH_WINDOW_LAG = SMOOTH_TIME_WINDOW/2;
const int CLEAR_SPATIAL_SMOOTH = 50;
const int COLLATED_SMOOTH_LAG = 3;

const float T_SMOOTH_COLLATED = 0.4; //0.2
const float T_TL0 = 0.3;
const float T_DIAG = 2*DELTA_SST*sqrt(FILTER_TIME_SIZE);
const float T_NN = 2*DELTA_SST;
const float MIN_TEMP = 271.15;
//const float MIN_TEMP = 270;
const float T_COLD_DT = -5;
const float T_WARM_DT = 10;
const float PASS2 = DELTA_SST;
const float T_INTERP = 0.5; //change to .3


const int COLLATED_LAG = 6;
const float GAMMA = COLLATED_LAG/2.0;

const float MU_CLEAR = 0.6;
const float MU_APPROX = 0.4;

const int LAND_KERNEL = 7;

const float T_DERIV = 0.04;
const float T_EIGEN = 0.1;
const float T_RATIO = 0.3;
const float T_COLD = 0.2;

const int SPACE_LAG = 9;
const int PASS_THRESH = 1000;
const int INTERP_DIST = 7;
