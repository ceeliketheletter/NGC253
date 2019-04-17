import os

#activate secondaries code

os.system("g++ `root-config --glibs` `root-config --cflags` -o secondaries_code.out secondaries_code.cpp")
os.system("./secondaries_code.out")
