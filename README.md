The code implements a method (constrained non-negative matrix factorization, CNMF) for simultaneous source extraction from large scale calcium imaging movies. The algorithm is presented in more detail here: http://arxiv.org/abs/1409.2903 
The CNMF code is written by Eftychios A. Pnevmatikakis and Liam Paninski, with inputs from Weijian Yang and Darcy S. Peterka. The demo script are written by Eftychios A. Pnevmatikakis, Liam Paninski, Weijian Yang and Darcy S. Peterka, 2015.

Use "demo_script2.m" to run the CNMF algorithm.
Use "dataProcessing.m" to extract results after "demo_script2.m" is ran. 
"dataProcessing.m" can save a .mat file which contains the CNMF parameters and the results.
Use "file_dataProcessing.m" to reload the results from the saved .mat file.
