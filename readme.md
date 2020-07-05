Created by Özgür Kara for EE242 Project
-----------------------------------------------------
You should compile the cpp file with writing into terminal "g++ -std=c++11 main.cpp Matrix.cpp -o main.exe". (use -std=c++11 because i used array initialization 
which is available in c++11) It creates an executable file. Then you can run this exe file, on the command prompt with writing 
"main.exe input.txt 1e-6 output.txt" where input.txt is the file that contains matrix and 1e-6 is the tolerance value and output.txt is the file that
eigenvalues and eigenvectors will be written into.

This program basically reads a square matrix from an input file and computes it's 2 dominant eigenvalue and corresponding eigenvector using
normalized power iteration algorithm and deflation method.

Program's approach to exceptions:
-----------------------------------------------------
If the input is zero matrix, then the output will be "Input is zero matrix!" written on the output.txt.

If the input is 1x1 matrix, then it has 1 eigenvalue and equals to its value. So there is 1 eigenvalue and eigenvector.

And because of the normalized iteration algorithm and deflation method, eigenvalues should be different from each other.
