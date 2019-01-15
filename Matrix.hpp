//
// Created by sin_een on 1/14/19.
//

#ifndef EX3_MATRIX_HPP
#define EX3_MATRIX_HPP

#include <iostream>
#include <vector>
#include "Complex.h"

using namespace std;

template <class T>
class Matrix
{
private:
    unsigned int columns;
    unsigned int rowNumbers;
    unsigned int size;
    vector<T> matrix;

public:

    /**
     * empty constructor that building a matrix of one
     */
    Matrix(): columns(1), rowNumbers(1), size(1)
    {
        matrix.push_back(T());
    }

    Matrix(unsigned int rows, unsigned int cols);

    Matrix(const Matrix<T>& other);

    Matrix(unsigned int rows, unsigned int cols, const vector<T>& cells);

    ~Matrix() {}

    /**
     * @param other
     * @return
     */
    Matrix<T>& operator=(const Matrix &other);

    /**
     * subtracts two matrixes  and puts the vlaues  on teh matrix in the class
     * @param other a second matrix
     * @return  matrix that is sum of both matrixes
     */
    Matrix<T> operator+(const Matrix &other);

    /**
     * multiplications two matrixes  and puts the vlaues  on teh matrix in the class
     * @param other a second matrix
     * @return matrix that is the subtraction between th e2 matrixes
     */
    Matrix<T> operator-(const Matrix &other);

    /**
     * addes the object of this class to other adn makes a new matrix and returns that
     * @param other another matrix3D
     */
    Matrix<T> operator*(const Matrix &other);

    /**
     * returns true if the matrix and the other have the same values
     */
    bool operator==(const Matrix& other);

    /**
     * returns true if the matrix and the other have the same values
     */
    bool operator!=(const Matrix& other);

    /**
     * only works for square matrixes ( colum = row) other wise retrun exception
     * @return a new matrix that isi trasposed to our matrix
     */
    Matrix<T> trans();


    /**
     * @return true if the matrix is square
     */
    bool isSquareMatrix();


    /**
     * prints out the matrix
     */
	template <class T0>
     friend ostream& operator<<(ostream& steam, const Matrix<T0> matrix1);

     /**
      *
      */
     T& operator()(unsigned int row, unsigned int col);

     /**
      *
      */
     T operator()(unsigned int row, unsigned int col) const;

typedef typename std::vector<T>::const_iterator const_iterator;
/// \returns the begin iterator
    const typename vector<T>::const_iterator begin()
    {
        return matrix.begin();
    }

/// \returns the end iterator
    const typename vector<T>::const_iterator end()
    {
        return matrix.end();
    }

    unsigned int rows()const;

    unsigned int cols()const;

    unsigned int getSize() const;

    vector<T> getMatrix() const;

    vector<T>& getMatrix();

    vector<T> getColumnIn(const Matrix<T>& m, unsigned int index);

    vector<T> getRowIn(const Matrix<T>& m, unsigned int index);

	static T multiplyColRow(vector<T> col, vector<T>row);
};

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols)
{
    this->size = rows*cols;
    this->rowNumbers = rows;
    this->columns = cols;
    this->matrix.clear();
    for(unsigned int i = 0; i < size; i++ )
    {
        matrix.push_back(T());
    }
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& other)
{
    this->rowNumbers = other.rows();
    this->columns = other.cols();
    this->size = rowNumbers * columns;

    this->matrix.clear();
    vector<T> m = other.getMatrix();
    for(unsigned int i = 0; i < size; i++ )
    {
        matrix.push_back(m.at(i));
    }
}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const vector<T> &cells)
{
    this->size = rows*cols;
    this->rowNumbers = rows;
    this->columns = cols;
    for(unsigned int i = 0; i < size; i++ )
    {
        matrix.push_back(cells.at(i));
    }
}


template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &other)
{
    this->size = other.getSize();
    this->rowNumbers = other.rows();
    this->columns = other.cols();
    this->matrix.clear();
    vector<T> m = other.getMatrix();
    for(unsigned int i = 0; i < size; i++ )
    {
        matrix.push_back(m.at(i));
    }
    return (*this);
}


template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix &other)
{
    unsigned int sizes = this->getSize();
    Matrix<T> sum = Matrix(this->rows(), this->cols());
    for(unsigned int i = 0; i < sizes; i++ )
    {
        sum.getMatrix().at(i) = this->getMatrix().at(i) + other.getMatrix().at(i);
    }
//	for( unsigned int r = 0; r < this->rows(); r ++ )
//	{
//		for( unsigned int c = 0; c < this->cols(); c ++ )
//		{
//			sum(r,c) = this(r, c) + other( r, c);
//		}
//	}


    return sum;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix &other)
{
    Matrix<T> sum = Matrix(this->rows(), this->cols());
    for(unsigned int i = 0; i < this->getSize() ; i++ )
    {
        sum.getMatrix().at(i) = this->getMatrix().at(i) - other.getMatrix().at(i);
    }
//	for( unsigned int r = 0; r < this->rows(); r ++ )
//	{
//		for( unsigned int c = 0; c < this->cols(); c ++ )
//		{
//			sum(r,c) = this(r, c) - other( r, c);
//		}
//	}

	return sum;
}

template<class T>
unsigned int Matrix<T>::rows()const
{
    return this->rowNumbers;
}

template<class T>
unsigned int Matrix<T>::cols()const
{
    return this->columns;
}

template<class T>
unsigned int Matrix<T>::getSize() const
{
    return this->size;
}

template<class T>
vector<T> Matrix<T>::getMatrix() const
{
    return this->matrix;
}

template<class T>
vector<T>& Matrix<T>::getMatrix()
{
    return this->matrix;
}


template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix &other)
{
	unsigned int newRows = this->rows();
	unsigned int newCols = other.cols();
	Matrix<T> mul = Matrix(newRows, newCols);
	for( unsigned int r = 0; r < newRows; r ++ )
	{
		vector<T> row = getRowIn(*this, r);
		for( unsigned int c = 0; c < newCols; c ++ )
		{
			vector<T> col = getColumnIn(other, c);
			mul(r,c) = multiplyColRow(col, row);
		}
	}
	return mul;
}

template<class T>
T &Matrix<T>::operator()(unsigned int row, unsigned int col)
{
    unsigned int index = (row * this->cols()) + col;
    return (*this).matrix.at(index);
}

template<class T>
T Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
    unsigned int index = (row * this->columns) + col;
    return this->matrix.at(index);
}

template<class T>
vector<T> Matrix<T>::getColumnIn(const Matrix &m,unsigned int index)
{
	if ( index < m.cols() and index >= 0)
	{
		vector<T> col;
		for (unsigned int i = 0; i < m.rows(); i ++)
		{
			col.push_back(m.matrix.at(index + (i * m.cols())));
		}
		return col;
	}
	else{
		throw "Out of bounderies"; //TODO exeption
	}
}

template<class T>
vector<T> Matrix<T>::getRowIn(const Matrix &m, unsigned int index)
{
	if (index < m.rows() and index >= 0 )
	{
		vector<T> col;
		for (unsigned int i = 0; i < m.cols(); i ++) {
			col.push_back(m.matrix.at((index * cols())+ i));
		}
		return col;
	}
	else{
		throw "Out of bounderies"; //TODO exeption
	}
}

template<class T>
T Matrix<T>::multiplyColRow(vector<T> col, vector<T> row)
{
	unsigned long sizeC = col.size();
	unsigned long sizeR = row.size();
	if( sizeC == sizeR)
	{
		T sum = T();
		for(unsigned long i = 0 ; i < sizeC; i++)
		{
			sum += col.at(i)* row.at(i);
		}
		return sum;
	}
	else
	{
		throw "colums and rows not equal cant multiply"; //TODO EXEPTIOM
	}
}
template<class T>
bool Matrix<T>::operator==(const Matrix &other)
{
	if (this->getSize() != other.getSize())
	{
		return false;
	}
	if (this->rows() != other.rows())
	{
		return false;
	}
	if (this->cols() != other.cols())
	{
		return false;
	}
	vector<T> vector1 = this->getMatrix();
	vector<T> vector2 = other.getMatrix();
	for (unsigned int i = 0 ; i < this->size; i++)
	{
		if( vector1.at(i) != vector2.at(i) )
		{
			return false;
		}
	}
	return true;
}

template<class T>
bool Matrix<T>::operator!=(const Matrix &other)
{
//	 return ! (this == *other);
	if (this->getSize() != other.getSize())
	{
		return true;
	}
	if (this->rows() != other.rows())
	{
		return true;
	}
	if (this->cols() != other.cols())
	{
		return true;
	}
	vector<T> vector1 = this->getMatrix();
	vector<T> vector2 = other.getMatrix();
	for (unsigned int i = 0 ; i < this->size; i++)
	{
		if( vector1.at(i) != vector2.at(i) )
		{
			return true;
		}
	}
	return false;
}

template<class T>
bool Matrix<T>::isSquareMatrix()
{
	return this->cols() == this->rows();
}

template<class T>
Matrix<T> Matrix<T>::trans()
{
	if (this->isSquareMatrix())
	{
		unsigned int length = this->rows();
		Matrix<T> transposed = Matrix(length, length);
		for( unsigned int r = 0; r < length; r ++ )
		{
			for (unsigned int c = 0; c < length; c ++)
			{
				transposed(r, c) = (*this)(c, r);
			}
		}
		return transposed;
	}
	else
	{
		//throw exception; //TODO exception
		throw "Out of bounderies";
	}
}
template <>
Matrix<Complex> Matrix<Complex> :: trans()
{
	if (this->isSquareMatrix())
	{
		unsigned int length = this->rows();
		Matrix<Complex> transposed = Matrix(length, length);
		for( unsigned int r = 0; r < length; r ++ )
		{
			for (unsigned int c = 0; c < length; c ++)
			{
				transposed(r, c) = (*this)(c, r).conj();

			}
		}
		return transposed;
	}
	else
	{
		throw "not square"; //TODO exception
	}
}

template<class T>
ostream &operator<<(ostream &stream, const Matrix<T> matrix1)
{
	unsigned long sizeC = matrix1.cols();
	unsigned long sizeR = matrix1.rows();
	for( unsigned int r = 0; r < sizeC; r ++ )
	{
		for (unsigned int c = 0; c < sizeR; c ++)
		{
			stream << matrix1(r,c) << "\t";
		}
		stream << endl;
	}
	return stream;
}




#endif //EX3_MATRIX_HPP
