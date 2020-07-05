#include <iostream>
#include "Matrix.h" //including Matrix class
#include <string>
#include <stdlib.h>
#include <fstream>

using namespace std;

//helper function that detects if matrix is zero matrix
bool is_zero_matrix(double** arr, int row)
{
	for(int i=0; i<row; i++)
	{
		for(int j=0;j<row;j++)
		{
			if(arr[i][j] != 0)
			{
				return 0;
			}
		}
	}
	return 1;
}

//function that performs deflation method and returns the reduced form (H A H-1)
Matrix deflation(int row, Matrix xk, Matrix A)
{
	Matrix I, v, H, similar;
	I = Matrix::Identity(row); //creates and identity matrix
	if(xk.get_num() >= 0) // if it's first element is positive, the sign is positive (to perform householder transformation)
	{
		v = xk + (xk.find_length_vector() * I.ith_col(0));
	}
	else // else the sign should be negative
	{
		v = xk - (xk.find_length_vector() * I.ith_col(0));
	}
	
	H = I - (2/(Matrix::multiply(v.transpose(),v).get_num())) * Matrix::multiply(v, v.transpose()); //performs householder transformation
	similar = Matrix::multiply(Matrix::multiply(H,A),H.transpose()); //finding similar matrix H A H-1 (H-1 = H.transpose())
	return similar.reduce(); //then reduces it to matrix and returns it (delete first row and first column)
	
}

int main(int argc, char** argv) 
{
	ofstream myfile; //creating osftream object
	Matrix A, x0, b, xk; //creating matrices
	double tolerance, ev1;
	int row,col = 0;
	string txt_file_A = argv[1]; //get first argument input name
	tolerance = atof(argv[2]); //get second argument tolerance
	A.read_A(txt_file_A); //set values of A matrix by reading txt file
	
	A.get_size(&row,&col); //get sizes of matrix
	myfile.open (argv[3]); //getting and opening output file
	
	if(is_zero_matrix(A.get_array(), row)) // if matrix is zero, no eigenvalues
	{
		myfile << "Input is zero matrix!";
		return 0;
	}
	
	if(row == 1) //if given matrix is 1x1, the eigenvalue is it's element
	{
		myfile << "Eigenvalue#1: "  <<  A.get_num() << endl;
		myfile << 1 << endl;
		return 0;
	}
	

	for(int k=0;k<2;k++) // it works for 2 step because we must find the dominant 2 eigenvalues.
	{
		//creating initial vector x0 [2,2,....]
		double** x = new double*[row-k]; 
		for(int i = 0; i < row-k; ++i)
		{
	    	x[i] = new double[1];
	    	x[i][0] = 2;
		}
		
		x0.set_values(row-k, 1, x); //x0 is initial vector	
		
		do{
			Matrix yk = Matrix::multiply(A,x0); // performs Ax0 = yk
			
			xk = (1/yk.find_max()) * yk; // performs yk / find_max(yk)   --> find_max(yk) converges to eigenvalue and xk converges to eigenvector
			ev1 = yk.find_max(); // it converges to eigenvalue

			x0 = xk; //then reassign the x0 parameter as xk
			Matrix diff = Matrix::multiply(A,xk) - ev1 * xk; //find the error vector
			
			if(diff.find_length_vector() < tolerance) //perform these operations until error < tolerance
			{
				myfile << "Eigenvalue#" << k+1 << ": " << ev1 << endl; // write the output eigenvalue
				if(k==0) // if it's first eigenvalue, write down the eigenvector to the corresponding this eigenvalue to the file
				{
					for(int o=0;o<row;o++)
					{
						myfile << xk.get_array()[o][0] << endl;	//write the eigenvector to the file			
					}
				}
				break;
			}
		}while(1); //in this while loop, it performs NORMALIZED POWER ITERATION until eigenvector and eigenvalue converges
		
		if(k==0) //perform DEFLATION operation 1 time.
		{
			A = deflation(row, xk, A);
		}
	}

	return 0;		
}

