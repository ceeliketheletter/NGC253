ofstream myfile;
myfile.open ("example.txt");  //This code creates a file called example.txt and inserts the for loop into it in the same way we are used to do with cout, but using the file stream myfile instead.
     
    
  
    
    
for (loggamma=.1; loggamma<=11; loggamma+=logbinsize)
  {
    ne_ss= NE_finite_time(& loggamma, & normalization);
    
    cout << loggamma << "\t" << ne_ss << endl; //to save to a file, use myfile instead of cout
  }
myfile.close(); //When we are finished with our input and output operations on a file we shall close it so that the operating system is notified and its resources become available again.
    
    
