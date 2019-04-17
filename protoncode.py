import os


#activate version 1
'''
os.system("g++ `root-config --glibs` `root-config --cflags` -o make_prot_dist_steady_state_w_args_version1.out make_prot_dist_steady_state_w_args_version1.cpp")
os.system("./make_prot_dist_steady_state_w_args_version1.out")
'''


#activate version 2: non normalized proton code
'''
os.system("g++ `root-config --glibs` `root-config --cflags` -o make_prot_dist_steady_state_w_args_version2.out make_prot_dist_steady_state_w_args_version2.cpp")
os.system("./make_prot_dist_steady_state_w_args_version2.out")
'''

#activate version 3: normalized proton code                                 

os.system("g++ `root-config --glibs` `root-config --cflags` -o normalized_proton_code.out normalized_proton_code.cpp")
os.system("./normalized_proton_code.out")





#activate the inputs to the original code

'''
os.system("g++ `root-config --glibs` `root-config --cflags` -o original_code.out original_code.cpp")
os.system("./original_code.out 400 2.0 2.1 0.0 0.0 0.0 15.0")
'''
