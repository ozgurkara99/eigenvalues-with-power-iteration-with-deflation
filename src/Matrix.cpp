#include <iostream> 
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"

using namespace std; 

//constructor function which has default parameters are row_num=1, col_num=1 as written in header file,
//when a Matrix object is created, it has dimensions (row_num and col_num) and double array that has 0 values inside of it. 
//Matrix a(n,m);  --> it creates a matrix object that has size (nxm) and all elements are 0.
//Matrix b;  --> it creates a matrix object that has size (1x1) and it's element is 0.
Matrix::Matrix (int row_num, int col_num)
{
	this->col_num = col_num;
	this->row_num = row_num;
    this->arr = new double*[row_num];
	for(int i = 0; i < row_num; ++i)
	    this->arr[i] = new double[col_num]();
	
}

//overloading summation operation of matrix. It sums the values of matrices elementwise.
//Matrix a, b;
//Matrix c = a + b; (summing all the values elementwise and returns a Matrix object)
Matrix operator+(const Matrix &x, const Matrix &y) 
{
    Matrix z(x.row_num, x.col_num);
    for(int i=0; i<x.row_num; i++)
    {
    	for(int j=0; j<x.col_num; j++)
    	{
    		z.arr[i][j] = x.arr[i][j] + y.arr[i][j];
		}
	}
	return z;
}

//overloading multiplication operation. You can multiply a matrix with scalar number
//double a; 
//Matrix b;
//Matrix c = a * b;  (it multiplies all the elements of b with a and returns another Matrix object) 
Matrix operator*(const double &x, const Matrix &y) 
{
    Matrix z(y.row_num, y.col_num);
    for(int i=0; i<y.row_num; i++)
    {
    	for(int j=0; j<y.col_num; j++)
    	{
    		z.arr[i][j] = x * y.arr[i][j];
		}
	}
	return z;
}

//overloading subtraction operation of matrix. It suubtracts the values of matrices elementwise.
//Matrix a, b;
//Matrix c = a - b; (subtracting all the values elementwise and returns a Matrix object)
Matrix operator-(const Matrix &x, const Matrix &y) 
{
    Matrix z(x.row_num, x.col_num);
    for(int i=0; i<x.row_num; i++)
    {
    	for(int j=0; j<x.col_num; j++)
    	{
    		z.arr[i][j] = x.arr[i][j] - y.arr[i][j];
		}
	}
	return z;
}


//function for reading a 2d array from A.txt and assign the values to Matrix object's values.
//function takes the filename in which matrix will be readed.
void Matrix::read_A(string filename)
{
	//open the file with the name "name" (converting it to const char array from string using c_str())
	ifstream file(filename.c_str());
	//initialize column number as zero
	int C = 0;
	int R = 0;
	//create null double 2d array (2d array is actually a pointer of pointer)
	double **arr = NULL;
	string line, temp;
	//while it gets line from the file, it updates the matrix 
	while (getline(file, line)) {
		//reallocate a pointer to pointer with the size of updating row number 
		arr = (double **)realloc(arr,(R+1)*sizeof(double*));
		//initialize it's column as null so that we can update (reallocate) when it finds a number in the line string
		arr[R] = NULL;
		//we should convert the string line to char array using c_str(). it converts to const char array so we duplicate it to make it char array
	    char * dup_str = strdup(line.c_str()); 
	    // Declaration of delimiter with any whitespace (tab or space)
	    //strtok splits the string using delimiters. it returns the first item until it sees the delimiter and the second part is in it's memory
		//for ex dup_str is 123 456 789, then tok is "123", and "456 789" is in the memory of strtok so that we can return other tokens in the while loop 
	    char* tok = strtok(dup_str, "\t "); 	  
	    while (1) {
	    	//if no token found, then break
	        if(tok == NULL)
	        	break;
	        //reallocate the columns of array when it finds a new number to add in a row
			arr[R] = (double*)realloc(arr[R],(C+1)*sizeof(double)); 
			//strtod converts the token to double and assign it to array
			arr[R][C] = strtod(tok,NULL);
			//update column number by incrementing
    		C = C + 1;
    		//find the next number
    		tok = strtok(NULL, "\t ");
	    } 
	    this->col_num = C;
   		C = 0;
		R = R + 1;		
  	}
  	this->row_num = R;
  	this->arr = arr;

}

//performs matrix multiplication (AX = B) and returns B matrix.
//Matrix A, X, B;
//B = Matrix::multiply(A,X);
Matrix Matrix::multiply(const Matrix& A, const Matrix& x)
{
	Matrix b(A.row_num, x.col_num);
	    
	for(int i = 0; i < A.row_num; ++i)
	{
		for(int j = 0; j < x.col_num; ++j)
		{
			for(int k=0; k<A.col_num; ++k)
			{
				b.arr[i][j] += A.arr[i][k] * x.arr[k][j];
			}
		}
	}
	return b;
}

//creates an Identity matrix with given size (n x n)
//Matrix I = Matrix::Identity(5)  --> will create an identity matrix with 5x5 size.
Matrix Matrix::Identity(int n)
{
	Matrix c(n,n);
	for(int i=0;i<n;i++)
	{
		c.arr[i][i] = 1;
	}
	return c;
}

//function for displaying the "row number", "column number" and "infinity norm" of a Matrix object
void Matrix::print()
{
	for(int i=0;i<row_num;i++)
	{
		for(int j=0;j<col_num;j++)
		{
			cout << arr[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "row size: " << row_num << "\n";
	cout << "col size: " << col_num << "\n";
	cout << "inf norm: " << this->find_inf_norm() << "\n";
}

//getter function to access row number and column number of a matrix object
void Matrix::get_size(int* row, int* col)
{
	*row = row_num;
	*col = col_num;	
}

//it returns the element that is in the ath row and bth column of matrix (starts from 0)
//default parameters are a=0, b=0
double Matrix::get_num(int a, int b)
{
	return arr[a][b];
}

//setter function for Matrix object. Sets the row number, column number and values manually if it's not read from any txt file or anything else.
void Matrix::set_values(int row_num, int col_num, double** arr)
{
	this->col_num = col_num;
	this->row_num = row_num;
	this->arr = arr;
}

//finding infinity norm of matrix. It returns the value.
//infinity norm is defined as max(absolute of sums of rows)
double Matrix::find_inf_norm()
{
	double inf = 0;
	double max = 0;
	double row_sum = 0;
	for(int i=0;i<row_num;i++)
	{
		row_sum = 0;
		for(int j=0;j<col_num;j++)
		{
			row_sum += arr[i][j];
		}
		row_sum = abs(row_sum);
		if(row_sum > max)
		{
			max = row_sum;
		}
		
	}
	return max;
}

//finding the sum of row value that value satisfies max(absolute of sum of rows)
//difference from inf norm is it's returning not the absolute value of the sum of rows value that has max magnitude.
//it's for finding the exact eigenvalue component.
double Matrix::find_max()
{
	double inf = 0;
	double max = 0;
	double row_sum = 0;
	for(int i=0;i<row_num;i++)
	{
		row_sum = 0;
		for(int j=0;j<col_num;j++)
		{
			row_sum += arr[i][j];
		}
		
		if(abs(row_sum)> abs(max))
		{
			max = row_sum;
		}
		
	}
	return max;
}

//function for finding the euclidian norm of a vector (nx1) size. it returns the value.
double Matrix::find_length_vector()
{
	double summ=0;
	for(int i=0;i<row_num;i++)
	{
		summ += arr[i][0] * arr[i][0];
	}
	return sqrt(summ);
}

//function for transpose operation, it returns a matrix object
//it returns another matrix object.
//Matrix a, b;
//b = a.transpose(); --> b is the transpose of a.
Matrix Matrix::transpose()
{
	Matrix c(this->col_num, this->row_num);
	for(int i=0; i<row_num;i++)
	{
		for(int j=0;j<col_num;j++)
		{
			c.arr[j][i] = this->arr[i][j];
		}
	}
	return c;
}

//function that returns a matrix object that is exactly the i'th column of given matrix.
//Matrix a,b;
//b = a.ith_col(3);  --> b is the vector that equals to 3th column of matrix a. (starts from 0)
Matrix Matrix::ith_col(int i)
{
	Matrix c(row_num, 1);
	for(int j=0;j<row_num;j++)
	{
		c.arr[j][0] = this->arr[j][i];
	}
	return c;
}

//function that returns (n-1 x n-1) size matrix from given matrix (n x n). It deletes the first row and first column 
//and returns the remaining part as a matrix object
//Matrix a,b;
//a = b.reduce(); --> it reduces b matrix and returns a matrix object that is assigned to a.
Matrix Matrix::reduce()
{
	Matrix c(row_num - 1 , col_num - 1);
	for(int i=0; i<c.row_num;i++)
	{
		for(int j=0; j<c.col_num;j++)
		{
			c.arr[i][j] = arr[i+1][j+1];
		}
	}
	return c;
}

//getter function that returns the values of matrix object
double** Matrix::get_array()
{
	return arr;
}

