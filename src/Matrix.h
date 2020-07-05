#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <string>

using namespace std;

//matrix class is created so that we can define matrix objects that have size variables, the value of array and matrix operations as functions.

class Matrix
{
	private:
		int col_num; //column number of matrix
		double** arr; //values of 2D matrix
		int row_num; //row number of matrix
		double find_inf_norm();	 //founding infinite norm
			
	public:
		static Matrix multiply(const Matrix& A, const Matrix& x);
		static Matrix Identity(int n);
		
		friend Matrix operator+(const Matrix &x, const Matrix &y);
		friend Matrix operator*(const double &x, const Matrix &y);
		friend Matrix operator-(const Matrix &x, const Matrix &y);
		
		double get_num(int a=0, int b=0);
		double find_length_vector();
		double find_max();
		
		double** get_array();		
				
		Matrix ith_col(int i=0);
		Matrix(int row_num=1, int col_num=1);
		Matrix transpose();
		Matrix reduce();
				
		void set_values(int row_num, int col_num, double** arr);
		void read_A(string name);
		void print();
		void get_size(int* row, int* col);
		


};

#endif


