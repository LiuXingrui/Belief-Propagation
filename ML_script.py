import os
#p_list=[0.0001,0.0004,0.0007,0.001,0.004,0.007,0.01,0.04,0.07,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.999]
#p_list=[0.01,0.04,0.07,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.999]
p_list=[0.1]
command= "rm data/ML_data.data "
os.system(command)
for p in p_list:
    command= "./MLSR ../S_BP/checkerboard_Hx_d3  ../S_BP/checkerboard_Hz_d3 %f 100 0 data/ML_data.data"%p
    os.system(command)
