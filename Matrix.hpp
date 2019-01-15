//
// Created by sin_een on 1/14/19.
//

#ifndef EX3_MATRIX_HPP
#define EX3_MATRIX_HPP

#define NOT_A_SQUARE "Matrix is not a Square matrix no trans"
#define NEXT_VALUE "\t"
#define NEXT_LINE "\n"
#define OUT_OF_BOUNDARIES  "index is out of boundaries"
#define UNCOMPATABLE_MATRIX_SIZE  "the sizes of the two Matrix do not match, can not perform task"
#define UNCOMPATABLE_VECTOR_SIZE  "the sizes of the two Vectors do not match, can not perform multiplication"

#include <exception>
#include <iostream>
#include <vector>
#include "Complex.h"

using namespace std;


class Exceptions : public exception
{
	private:
		string message_Excep;

	public:
		Exceptions(const string& message): message_Excep(message){}

		~Exceptions() throw() {}

		virtual constexpr decltype(auto) what() const noexcept;

};

/**
 * @return the exeption msg
 */
constexpr decltype(auto) Exceptions::what() const
{
	return this->message_Excep;
}


/**
 *
 * @tparam T
 */
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
     * empty constructor that building a matrix of one cell that has zero in its value
     */
    Matrix(): columns(1), rowNumbers(1), size(1)
    {
        matrix.push_back(T());
    }

	/**
	 * constructor that building a matrix of size rows * colums note the vector / actual matrix is  empty at this level
	 * @param rows number of rows we want in our matrix
	 * @param cols number of columns in out matrix
	 */
    Matrix(unsigned int rows, unsigned int cols);

	/**
	 * a copy constructor that building a matrix with same values of the other matirx
	 * @param other  a matrix we want to copy
	 * note this deep copies and the new matrix is not related to the second one after the copy
	 */
    Matrix(const Matrix<T>& other);

	/**
	 * constructor that building a matrix of size rows * colums , and fills the matrix we have with the values in cells
	 * @param rows number of rows we want in our matrix
	 * @param cols number of columns in out matrix
	 * @param cells a vector with the values we want in the matrix
	 * throws an exception in case the number of values in cells are not equal to the number of cells in the matrix
	 * we want to build aka rows * columns
	 */
    Matrix(unsigned int rows, unsigned int cols, const vector<T>& cells);

	/**
	 * destructor = defulte since we do not allocate anything
	 */
    ~Matrix() {}

    /**
     * @param other puts the values of other in this matrix
     * @return refrence to this
     * throws an exception in case the sizes are not equal  aka rows * columns
     */
    Matrix<T>& operator=(const Matrix &other);

    /**
     * sums up the two matrixes  and puts the values in the matrix in the class
     * @param other a second matrix
     * @return  matrix that is sum of both matrixes
     * throws an exception in case the sizes are not equal  aka rows * columns
     *
     */
    Matrix<T> operator+(const Matrix &other);

	/**
	  * subtracts two matrixes  and puts the vlaues  on teh matrix in the class
	  * @param other a second matrix
	  * @return  matrix that is subtracts of both matrixes
	  * throws an exception in case the sizes are not equal  aka rows * columns
	  *
	  */
    Matrix<T> operator-(const Matrix &other);

	/**
	* multiplications two matrixes  and puts the vlaues  on teh matrix in the class
	* @param other a second matrix
	* @return matrix that is the subtraction between the 2 matrixes
	* catches an exception in case the the size of the rows of this is not equal to size of colums in other
	*/
    Matrix<T> operator*(const Matrix &other);

    /**
     * returns true if the matrix and the other have the same values
     */
    bool operator==(const Matrix& other);

    /**
     * returns true if the matrix and the other do not have the same values
     */
    bool operator!=(const Matrix& other);

    /**
     * only works for square matrixes ( colum = row) other wise throws exception
     * @return a new matrix that is the trasposed matrix to our matrix
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
      * returns a refrence of the cell in row col
      * throws exception if it was out of bounderes
      */
     T& operator()(unsigned int row, unsigned int col);

	/**
     * returns the cell of the matrix in (row ,col)  (look not change thus its const and doesnt return refrence)
     * throws exception if it was out of bounderes
     */
	T operator()(unsigned int row, unsigned int col) const;

typedef typename std::vector<T>::const_iterator const_iterator;
	/**

	 * @return the begin iterator
	 */
    const typename vector<T>::const_iterator begin()
    {
        return matrix.begin();
    }

	/**
	 *
	 * @return the end iterator
	 */
    const typename vector<T>::const_iterator end()
    {
        return matrix.end();
    }

	/**
	 * @return the number of rows in our matrix
	 */
    unsigned int rows()const;

	/**
 	* @return the number of columns in our matrix
 	*/
    unsigned int cols()const;

	/**
 	* @return the number of cells in our matrix
 	*/
    unsigned int getSize() const;

	/**
	 * @return gets the actual matrix we have Vector witht the values of the cells
	 */
    vector<T> getMatrix() const;

	/**
	* @return a refrnece of the actual matrix we have aka vector with the values in case we want to change
	 * values in the matrix
	*/
    vector<T>& getMatrix();

	/**
	 * @param m  the matrix we want the colum from
	 * @param index number of colum we wanted
	 * @return a vector that is that colum form the matrix m
	 * throws exeption if index is out of boundreis
	 */
    vector<T> getColumnIn(const Matrix<T>& m, unsigned int index);

	/**
 	 * @param m  the matrix we want the colum from
 	 * @param index number of row we wanted
 	 * @return a vector that is that row form the matrix m
 	 * throws exeption if index is out of boundreis
 	 */
    vector<T> getRowIn(const Matrix<T>& m, unsigned int index);

	/**
	 * @param col a colum of teh matrix
	 * @param row a row of teh matrix
	 * @return returns the value of the multiplication of two vectors
	 * throws exception in case the the size of the rows of this is not equal to size of colums in other
	 */
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
	if ( rows * cols == cells.size())
	{
		this->size = rows * cols;
		this->rowNumbers = rows;
		this->columns = cols;
		for (unsigned int i = 0; i < size; i ++)
		{
			matrix.push_back(cells.at(i));
		}
	}
	else
	{
		throw Exceptions(UNCOMPATABLE_MATRIX_SIZE);
	}
}


template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &other)
{
	if(other.getSize() != this->getSize())
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
	else
	{
		throw Exceptions(UNCOMPATABLE_MATRIX_SIZE);
	}
}


template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix &other)
{
	if(other.getSize() != this->getSize())
	{
		unsigned int sizes = this->getSize();
		Matrix<T> sum = Matrix(this->rows(), this->cols());
		vector<T> sumMatrix = sum.getMatrix();
		vector<T> thisMatrix = this->getMatrix();
		vector<T> otherMatrix = other.getMatrix();
		for (unsigned int i = 0; i < sizes; i ++)
		{
			sumMatrix.at(i) = thisMatrix.at(i) + otherMatrix.at(i);

		}
		return sum;
	}
	else
	{
		throw Exceptions(UNCOMPATABLE_MATRIX_SIZE);
	}
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix &other)
{
	if(other.getSize() != this->getSize())
	{
		Matrix<T> sum = Matrix(this->rows(), this->cols());
		vector<T> sumMatrix = sum.getMatrix();
		vector<T> thisMatrix = this->getMatrix();
		vector<T> otherMatrix = other.getMatrix();
		for(unsigned int i = 0; i < this->getSize() ; i++ )
		{
			sumMatrix.at(i) = thisMatrix.at(i) - otherMatrix.at(i);
		}

		return sum;
	}
	else
	{
		throw Exceptions(UNCOMPATABLE_MATRIX_SIZE);
	}
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
    return (*this).matrix;
}


template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix &other)
{
	try
	{
		unsigned int newRows = this->rows();
		unsigned int newCols = other.cols();
		Matrix<T> mul = Matrix(newRows, newCols);
		for (unsigned int r = 0; r < newRows; r ++)
		{
			vector<T> row = getRowIn(*this, r);
			for (unsigned int c = 0; c < newCols; c ++)
			{
				vector<T> col = getColumnIn(other, c);
				mul(r, c) = multiplyColRow(col, row);
			}
		}
		return mul;
	}
	catch (Exceptions msg)
	{
		throw;
	}
}

template<class T>
T &Matrix<T>::operator()(unsigned int row, unsigned int col)
{
	if (this->rows() > row & this->cols() > col )
	{
		unsigned int index = (row * this->cols()) + col;
		return (*this).matrix.at(index);
	}
	else if (this->rows() < row)
	{
		throw Exceptions(OUT_OF_BOUNDARIES);
	}
	else
	{
		throw Exceptions(OUT_OF_BOUNDARIES);
	}
}

template<class T>
T Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
	if (this->rows() > row & this->cols() > col )
	{
		unsigned int index = (row * this->cols()) + col;
		return this->matrix.at(index);
	}
	else if (this->rows() < row)
	{
		throw Exceptions(OUT_OF_BOUNDARIES);
	}
	else
	{
		throw Exceptions(OUT_OF_BOUNDARIES);
	}
}

template<class T>
vector<T> Matrix<T>::getColumnIn(const Matrix &m, unsigned int index)
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
	else
	{
		throw Exceptions(OUT_OF_BOUNDARIES);
	}
}

template<class T>
vector<T> Matrix<T>::getRowIn(const Matrix &m, unsigned int index)
{
	if (index < m.rows() and index >= 0 )
	{
		vector<T> col;
		for (unsigned int i = 0; i < m.cols(); i ++)
		{
			col.push_back(m.matrix.at((index * cols()) + i));
		}
		return col;
	}
	else
	{
		throw Exceptions(OUT_OF_BOUNDARIES);
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
		throw Exceptions(UNCOMPATABLE_VECTOR_SIZE);
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
		throw Exceptions(NOT_A_SQUARE);
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
		throw Exceptions(NOT_A_SQUARE);
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
			stream << matrix1(r, c) << NEXT_VALUE;
		}
		stream << NEXT_LINE;
	}
	return stream;
}




#endif //EX3_MATRIX_HPP
