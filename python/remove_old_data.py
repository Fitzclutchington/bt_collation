import os

data_folder = "../data"
granule_folders = [ "approx", "approx2","clear","collated_mat","collated_mat2","pass2","reinstated",
					"smooth","smooth_collate","TL0","TL1","TL2"]

for folder in granule_folders:
	current_folder = os.path.join(data_folder,folder)
	print current_folder
	for f in os.listdir(current_folder):
		os.remove(os.path.join(current_folder,f))
