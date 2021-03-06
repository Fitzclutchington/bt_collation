import os
import sys

if len(sys.argv) < 3:
    print "usage: python make_directory.py <output_folder> <brightness_temp>"
    print "brightness temp = 0 run for generated data"
    print "brightness temp = 1 run for brightness temps"
    sys.exit()

save_date = sys.argv[1]
data_folder = "../data"
top_folder = os.path.join(data_folder,save_date)

min_ys = [4600, 4150, 3500, 3500, 1300, 500]
max_ys = [5101, 4651, 4501, 4001, 1801, 1001]
min_xs = [2700, 3250, 1200, 3000, 1200, 2300]
max_xs = [3701, 4251, 1701, 4001, 2201, 3301]



if sys.argv[2] == '0':
    data = ['original_sst', 'acspo','pass2','reinstated','collated_mat2','smooth_collate']
else:
    data = ['brightness_temperature_08um6', 'brightness_temperature_10um4','brightness_temperature_11um2','brightness_temperature_12um3']

data_folder = os.path.join("../data",top_folder)

if not os.path.exists(top_folder):
    os.makedirs(top_folder)


for j in range(len(min_xs)):
    max_y = max_ys[j]
    min_y = min_ys[j]
    max_x = max_xs[j]
    min_x = min_xs[j]

    crop_folder = '_'.join(map(str,[min_y,max_y,min_x,max_x]))
    print crop_folder

    folder2 = '/'.join([top_folder,crop_folder])
    if not os.path.exists(folder2):
        os.makedirs(folder2)

    
    for d in data:
        folder3 = '/'.join([folder2,d])
        if not os.path.exists(folder3):
            os.makedirs(folder3)
    