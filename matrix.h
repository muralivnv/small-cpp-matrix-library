#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdarg>
#include <stdlib.h>

using mat2d = std::vector<std::vector<double> > ;

class Matrix{
		int row, column;
        mat2d M;
    public:
        Matrix(int a = 2, int b = 2):row(a),column(b){
            M.resize(row,std::vector<double>(column,0));}
		
		void val(double,...);       				   // For Assigning values into Matrix,
		
		Matrix operator*(const Matrix&);               // Matrix * Matrix
		friend Matrix operator*(double,const Matrix&); // scalar * Matrix
		friend Matrix operator*(const Matrix&,double); // Matrix * scalar
		
		Matrix t();                             	   // For doing Transpose
		double det();                                  // For doing Determinant
		void lu(Matrix&,Matrix&);	                   // For LU Decomposition
		Matrix inv();  								   // For doing Matrix inverse
        
		friend Matrix diag(int,...);                   // For creating diag_matrix,
		friend Matrix zeros(int, int);                 // For creating Zeros(rows,columns)
		friend Matrix ones (int, int);                 // For creating Ones(rows, columns) 
		friend Matrix eye(int, int);                   // For creating I(rows,columns)
		friend Matrix reshape(const Matrix&,int,int);  // For reshaping the matrix
		friend Matrix rand(int, int); 				   // For random matrix

		Matrix operator+(const Matrix&);               // For doing Matrix + Matrix
		friend Matrix operator+(double, const Matrix&);// For doing Scalar + Matrix
		friend Matrix operator+(const Matrix&,double); // For doing Matrix + Scalar
        
		Matrix operator-(const Matrix&);               // For doing Matrix - Matrix
		friend Matrix operator-(double, const Matrix&);// For doing Scalar - Matrix		
		friend Matrix operator-(const Matrix&,double); // For doing Matrix - Scalar
		
		Matrix& operator|(const Matrix&);			   // For horz concatenation [A;B]
		Matrix& operator||(const Matrix&);             // For vert concatenation [A B]
		
		int nrows();                                    // For getting no.of rows
		
		int ncols();                                    // For getting no.of column
		
		int nelem();                                   // For getting no.of elements
		std::vector<double> get_diag();                // For getting diag elements
		
		double operator()(int,int);                    // To Access individual elements
		void set(int,int,const double&);               // To set individual elements to a specific value
		friend std::ostream& operator<<(std::ostream& os, const Matrix&);  // For Printing

};


void Matrix::val(double a,...)
{
	va_list vals;
	va_start(vals,a);
	for(int i=0; i< row; i++)
	{
		for(int j=0; j<column; j++)
		{
		    if (i == 0 && j==0)
		    {
		        M[i][j] = a;
		    }
		    else 
			{
				M[i][j] = va_arg(vals,double);
		    }
		}
	}
}


Matrix Matrix::operator*(const Matrix& A)
{
  if (this->column != A.row)
	{
	  std::cerr<<"Dimensions inconsistent for Matrix Multiplication";
	  exit(1);
	}
  Matrix M2(this->row,A.column);
    
  for (int i=0;i<this->row;i++)
    {
	  for (int j=0;j<A.column;j++)
        {
		  for (int k=0;k<A.row;k++)
            {
			  M2.M[i][j]+=this->M[i][k]*A.M[k][j];
            }
        }
    }
  return M2;
    
}

Matrix operator*(double a,const Matrix& A)
{
	Matrix M2(A.row,A.column);
	for (int i=0; i<M2.row; i++)
	{
		for (int j=0; j<M2.column; j++)
		{
			M2.M[i][j] = a*A.M[i][j];
		}
	}
	return M2;
}

Matrix operator*(const Matrix& A,double a)
{
	return a*A;
}

Matrix Matrix::t()
{
    Matrix M2(column,row);    
    for(int i=0; i<M2.row; i++)
    {
        for(int j=0; j<M2.column; j++)
        {
            M2.M[i][j] = M[j][i];
        }
    }
    return M2;
}


double Matrix::det()
{
	if (column != row)
	{
		std::cerr<<"Determinant cannot be evaluted for Non-Squared Matrix";
		exit(1);
	}
    Matrix L(row,column);
    Matrix U(row,column);
    lu(L,U);
	double sumL = 1;
	double sumU = 1;
	for (int i=0;i<row;i++)
	{
			sumL = sumL*L.M[i][i];
			sumU = sumU*U.M[i][i];
	}
	// For erasing L and U memory,
	for (int i=0;i<L.row;i++)
	{
			L.M[i].erase(L.M[i].begin(),L.M[i].end());
			U.M[i].erase(U.M[i].begin(),U.M[i].end());
	}
	L.M.erase(L.M.begin(),L.M.end());
	U.M.erase(U.M.begin(),U.M.end());
	return (sumL*sumU);
	//
	
}

void Matrix::lu(Matrix& L,Matrix& U)
{
	if (column != row)
	{
		std::cerr<<"LU decomposition cannot be evaluted for Non-Squared Matrix";
		exit(1);
	}
    int i, j, k;
	double sum = 0;
    int n = row;
	for (i = 0; i < n; i++)
	{
	    L.M[i][0] = M[i][0];
		U.M[i][i] = 1;
	}

	for (j = 0; j < n; j++) 
	{
		for (i = j; i < n; i++) 
		{
			sum = 0;
			for (k = 0; k < j; k++) 
			{
				sum = sum + L.M[i][k] * U.M[k][j];	
			}
			L.M[i][j] = M[i][j] - sum;
		}

		for (i = j; i < n; i++) 
		{
			sum = 0;
			for(k = 0; k < j; k++) 
			{
				sum = sum + L.M[j][k] * U.M[k][i];
			}
			U.M[j][i] = (M[j][i] - sum) / L.M[j][j];
		}
	}
}

Matrix Matrix::inv()
{
	if (column != row)
	{
		std::cerr<<"Inverse cannot be evaluted for Non-Squared Matrix use pinv()";
		exit(1);
	}
    Matrix L(row,column);
    Matrix U(row,column);
    lu(L,U);
	
	mat2d I(row,std::vector<double>(column,0));
	mat2d D(row,std::vector<double>(column,0));
    for(int i=0;i<row;i++)
    {
    	I[i][i]=1;
    }
    
    double tmp = 0;
    for (int i=0;i<row;i++)
    {
        D[i][i] = 1/L.M[i][i];
    }
    // doing L*D = B->(I);
    for(int i = 0;i<column;i++)
    {
        for(int j=1;j<row;j++)
        {
            for(int k=0;k<=j-1;k++)
            {
                tmp+= L.M[j][k]*D[k][i]; 
            }
            D[j][i] = (I[j][i]- tmp)/L.M[j][j];
            tmp = 0;
        }
    }
    Matrix A_I(row,column);
    for(int i=0; i<column; i++)
    {
        A_I.M[row-1][i] = D[row-1][i];
    }

	for(int i=0;i<column;i++)
	{
	    for(int j=column-2;j>=0;j--)
	    {
		    for(int k=column-1; k>j; k--)
		    {
		    	tmp+= U.M[j][k]*A_I.M[k][i];
		    }
		    A_I.M[j][i] = D[j][i]-tmp;
		    tmp=0;
	    }
	}
	// For erasing L, U, I, D memory,
	/* for (int i=0;i<L.row;i++)
	{
			L.M[i].erase(L.M[i].begin(),L.M[i].end());
			U.M[i].erase(U.M[i].begin(),U.M[i].end());
			I[i].erase(I[i].begin(),I[i].end());
			D[i].erase(D[i].begin(),D[i].end());
	}
	L.M.erase(L.M.begin(),L.M.end());
	U.M.erase(U.M.begin(),U.M.end());
	I.erase(I.begin(),I.end());
	D.erase(D.begin(),D.end()); */
	return A_I;
}

Matrix diag(int num,...)
{
	Matrix M1(num,num);
    va_list num_list;
	va_start(num_list,num);
    for(int i = 0; i < num; i++)
    {
        M1.M[i][i] = va_arg(num_list,double);
    }
	va_end(num_list);
    return M1;
}


Matrix zeros(int a,int b)
{
	Matrix M1(a,b);
    return M1;
}

Matrix ones(int a,int b)
{
	Matrix M1(a,b);
	for(int i=0;i<M1.row;i++)
	{
		for(int j=0;j<M1.column;j++)
		{
			M1.M[i][j] = 1;
		}
	}
	return M1;
}

Matrix eye(int a,int b)
{
   Matrix M1(a,b);
    for(int i=0;i<a;i++)
    {
    	for(int j=0;j<b;j++)
    	{
    		if(i==j)
    		{
    			M1.M[i][j]=1;
    		}
    		else
    		{
    			M1.M[i][j]=0;
    		}
    	}
    }
    return M1;
}

Matrix reshape(const Matrix& A, int a,int b)
{
	if (A.row*A.column != a*b)
	{
		std::cerr<<"check the reshaping dimensions-> dimensions inconsistent";
		exit(1);
	}
	Matrix M1(a,b);
	int n_elem = a*b; //no.of elements,
	mat2d tmp(n_elem,std::vector<double> (1,0));

	int n=0;
	for(int i=0; i<A.column; i++)
	{
		for(int j=0; j<A.row; j++)
		{
			tmp[n][0] = A.M[j][i];
			n+=1;
		}
	}
	n=0;
// Assigning the values to new matrix,
		for(int i=0; i<b; i++) // columns
		{
			for(int j=0; j<a; j++) //rows
			{
				M1.M[j][i] = tmp[n][0];
				n+=1;
			}
		}
		for (unsigned int i=0;i<tmp.size();i++)
		{
			tmp[i].erase(tmp[i].begin(),tmp[i].end());
		}
	tmp.erase(tmp.begin(),tmp.end());
	return M1;
}


Matrix rand(int a, int b)
{
	Matrix M1(a,b);
    srand(time(NULL));
    for(int i=0;i<M1.row;i++)
    {
        for(int j=0;j<M1.column;j++)
        {
            M1.M[i][j]= rand() % 10 + 1;
			M1.M[i][j]*= 0.1;
        }
    }

    return M1;
}

Matrix Matrix::operator+(const Matrix& A)
{
	if (this->column != A.column || this->row != A.row)
	{
		std::cerr<<"Matrices cannot be added, Dimensions mismatch";
		exit(1);
	}
	Matrix M1(this->row,A.column);    
    for (int i=0;i<this->row;i++)
    {
        for (int j=0;j<A.column;j++)
        {
                M1.M[i][j]=this->M[i][j]+A.M[i][j];    
        }
    }
    return M1;
}

Matrix operator+(double a, const Matrix& A)
{
	Matrix M1(A.row,A.column);
	
	for(int i=0; i<A.row;i++)
	{
		for(int j=0; j<A.column;j++)
		{
			M1.M[i][j] = a + A.M[i][j];
		}
	}
	return M1;
}

Matrix operator+(const Matrix& A,double a)
{
	return (a+A);
}

Matrix Matrix::operator-(const Matrix& A)
{
	if (this->column != A.column || this->row != A.row)
	{
		std::cerr<<"Matrices cannot be subtracted, Dimensions mismatch";
		exit(1);
	}
	Matrix M1(this->row,A.column);
    
    for (int i=0;i<this->row;i++)
    {
        for (int j=0;j<A.column;j++)
        {
                M1.M[i][j]=this->M[i][j]-A.M[i][j];    
        }
    }
    return M1;
}

Matrix operator-(double a, const Matrix& A)
{
	Matrix M1(A.row,A.column);
	
	for(int i=0; i<A.row;i++)
	{
		for(int j=0; j<A.column;j++)
		{
			M1.M[i][j] = a - A.M[i][j];
		}
	}
	return M1;
}

Matrix operator-(const Matrix& A,double a)
{
	return -a+A ;
}

Matrix& Matrix::operator|(const Matrix& C)
{
	if (this->column != C.column)
	{
		std::cerr<<"horz cat failed, Dimensions mismatch";
		exit(1);
	}
	int v_r = this->row;
    for(int i=0;i<C.row;i++)
    {
        this->M.push_back(std::vector<double> (C.column));
    }
    
    this->row+=C.row;
    
    for(int i=0;i<C.row;i++)
    {
        for(int j=0;j<C.column;j++)
        {
            this->M[i+v_r][j] = C.M[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator||(const Matrix& A)
{
	if (this->row != A.row)
	{
		std::cerr<<"Vertical cat failed, Dimensions mismatch";
		exit(1);
	}
    int v_c = this->column;
    for(int i=0;i<A.row;i++)
    {
        this->M[i].resize(this->column+A.column);
    }
    
    this->column+=A.column;
    
    for(int i=0;i<A.row;i++)
    {
        for(int j=0;j<A.column;j++)
        {
            this->M[i][j+v_c] = A.M[i][j];
        }
    }
    return *this;
}

int Matrix::nrows()
{ return row; }

int Matrix::ncols()
{ return column; }

int Matrix::nelem()
{ return row*column; }

std::vector<double> Matrix::get_diag()
{
	std::vector<double> d_elem(row,0);
	for(int i=0;i<row;i++)
	{
		d_elem[i] = M[i][i];
	}
	return d_elem;
}

double Matrix::operator()(int i,int j)
{
	return this->M[i][j];
}

void Matrix::set(int a,int b,const double& value)
{
	M[a][b] = value;
}

std::ostream& operator<<(std::ostream& os, const Matrix& A)  
{  
    for (int i=0;i<A.row;i++)
    {
        for(int j=0;j<A.column;j++)
        {    
            os << std::setw(12)<<A.M[i][j];
        }
        os<<"\n";
    }
    return os;  
} 

#endif

/* #define n 1000
#define BlockSize  100
int main()
{
	int a[n][n],b[n][n],c[n][n];
	c[0][0]=0;
	for( i1=0;i1<(n/BlockSize);++i1)
	{
		for(j1=0;j1<(n/BlockSize);++j1)
		{
			for(k1=0;k1<(n/BlockSize);++k1)
			{
				for(i=i1=0;i<min(i1+BlockSize-1);++i)
				{
					for(j=j1=0;j<min(j1+BlockSize-1);++j)
					{
						for(k=k1;k<min(k1+BlockSize-1);++k)
						{               
							c[i][j] = c[i][j] + a[i][k] * b[k][j]
						}
					}
				}
			}
		}
    }
 return 0;
} */

