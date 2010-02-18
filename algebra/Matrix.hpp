#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <vector>
#include <iosfwd>

template<class F>
class Matrix{
public:
    typedef F Base;
    Matrix<F>() 
    : mNumRows(0), mNumCols(0)
        {}
    
    Matrix<F>(size_t numRows, size_t numCols) 
    : mNumRows(numRows), mNumCols(numCols)
        {
            initVec();
        }
    
    F& operator () (size_t i, size_t j) { return mData[i][j]; }
    const F& operator () (size_t i, size_t j) const { return mData[i][j]; }
    
    bool
    operator < (const Matrix<F>& other) const{
	for (size_t i = 0; i < rows(); ++i)
	    for (size_t j = 0; j < cols(); ++j)
		if (mData[i][j] < other(i, j))
		    return true;
		else if (other(i, j) < mData[i][j])
		    return false;
	return false;
    }
    
    class Col{
    public:
        Col(size_t index, const Matrix& mat)
            : mLen(mat.mNumRows), mIndex(index), mMat(mat)
            {}
        
        const F& operator () (size_t i) const { return mMat(i, mIndex); }
        const size_t mLen;
        
    private:
        const size_t mIndex;
        const Matrix<F>& mMat;
    };
    
    class Row{
    public:
        Row(size_t index, const Matrix& mat)
            : mLen(mat.mNumCols), mIndex(index), mMat(mat)
            {}
        
        const F& operator () (size_t j) const { return mMat(mIndex, j); }
        const size_t mLen;
        
        F operator*(const Col& col) const
            {
                assert(mLen == col.mLen);
                F sum(0);
                const Row& row = *this;
                for (size_t i = 0; i < row.mLen; ++i)
                    sum = sum + row(i) * col(i);
                return sum;
            }
        
    private:
        const size_t mIndex;
        const Matrix<F>& mMat;
    };


    Col col(size_t i) const { return Col(i, *this); }
    Row row(size_t i) const { return Row(i, *this); }

    size_t rows() const { return mNumRows; }
    size_t cols() const { return mNumCols; }


    static Matrix<F> eye(size_t n)
        {
            Matrix<F> e(n, n);
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    e(i, j) = (i == j) ? F(1) : F(0);
            return e;
        }
    static Matrix<F> zeros(size_t n1, size_t n2)
        {
            Matrix<F> z(n1, n2);
            for (size_t i = 0; i < n1; ++i)
                for (size_t j = 0; j < n2; ++j)
                    z(i, j) = F(0);
            return z;
        }
    static Matrix<F> ones(size_t n1, size_t n2)
        {
            Matrix<F> z(n1, n2);
            for (size_t i = 0; i < n1; ++i)
                for (size_t j = 0; j < n2; ++j)
                    z(i, j) = F(1);
            return z;
        }


private:
    void initVec()
        {
            mData.resize(mNumRows);
            for (size_t i = 0; i < mNumRows; ++i)
                mData[i].resize(mNumCols);
        }


private:
    size_t mNumRows;
    size_t mNumCols;  
    std::vector<std::vector<F> > mData;
};

template<class F>
Matrix<F>
tposed(const Matrix<F>& m)
{
    Matrix<F> mt(m.cols(), m.rows());
    for (size_t i = 0; i < mt.rows(); ++i)
        for (size_t j = 0; j < mt.cols(); ++j)
            mt(i, j) = m(j, i);
    return mt;
}

template<class F>
Matrix<F>
operator * (const Matrix<F>& m1, const Matrix<F>& m2)
{
    assert(m1.cols() == m2.rows());
    Matrix<F> m(m1.rows(), m2.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            m(i, j) = m1.row(i) * m2.col(j);
    return m;
}

template<class F>
Matrix<F>
operator + (const Matrix<F>& m1, const Matrix<F>& m2)
{
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    Matrix<F> m(m1.rows(), m1.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            m(i, j) = m1(i, j) + m2(i, j);
    return m;
}

template<class F>
Matrix<F>
operator - (const Matrix<F>& m1, const Matrix<F>& m2)
{
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    Matrix<F> m(m1.rows(), m1.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            m(i, j) = m1(i, j) - m2(i, j);
    return m;
}

template<class F>
Matrix<F>
operator * (const Matrix<F>& m, const F& f)
{
    Matrix<F> r(m.rows(), m.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            r(i, j) = m(i, j) * f;
    return r;
}

template<class F>
Matrix<F>
operator % (const Matrix<F>& m, size_t p)
{
    Matrix<F> r(m.rows(), m.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            r(i, j) = ((m(i, j) % p) + p) % p;
    return r;
}

template<class F>
Matrix<F>
operator / (const Matrix<F>& m, const F& f)
{
    Matrix<F> r(m.rows(), m.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            r(i, j) = m(i, j) / f;
    return r;
}

template<class F>
Matrix<F>
operator * (const F& f, const Matrix<F>& m)
{
    Matrix<F> r(m.rows(), m.cols());
    for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
            r(i, j) = f * m(i, j);
    return r;
}

template<class F>
bool
operator == (const Matrix<F>& m1, const Matrix<F>& m2)
{
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    
    for (size_t i = 0; i < m1.rows(); ++i)
        for (size_t j = 0; j < m1.cols(); ++j)
            if (m1(i, j) != m2(i, j))
                return false;
    return true;
}

template<class F>
std::ostream& 
operator << (std::ostream& os, const Matrix<F>& m)
{
    for (size_t i = 0; i < m.rows(); ++i){
        for (size_t j = 0; j < m.cols(); ++j)
            os << std::setw(10) << m(i, j) << " ";
        os << std::endl;
    }
    return os;
}

template<class F>
void
divRow(Matrix<F>& m, size_t i, F x)
{
    for (size_t j = 0; j < m.cols(); ++j)
        m(i, j) = m(i, j) / x;
}

template<class F>
void
swapRows(Matrix<F>& m, size_t i1, size_t i2)
{
    for (size_t j = 0; j < m.cols(); ++j){
        F temp = m(i1, j);
        m(i1, j) = m(i2, j);
        m(i2, j) = temp;
    }
}

template<class F>
void
subRows(Matrix<F>& m, size_t i1, size_t i2, F x)
{
    for (size_t j = 0; j < m.cols(); ++j)
        m(i1, j) = m(i1, j) - m(i2, j) * x;
}

template<class F>
Matrix<F>
pinv(const Matrix<F>& m)
{
    Matrix<F> md(m.rows(), 2 * m.cols());
    for (size_t i = 0; i < md.rows(); ++i)
        for (size_t j = 0; j < md.cols() / 2; ++j)
            md(i, j) = m(i, j);
    for (size_t i = 0; i < md.rows(); ++i)
        for (size_t j = md.cols() / 2; j < md.cols(); ++j)
            md(i, j) = 
                (i == j - md.cols() / 2) ? F(1) : F(0);
    for (size_t iter = 0; iter < min(m.rows(), m.cols()); ++iter){
        size_t ipiv = iter;
        for (; ipiv < md.rows() && md(ipiv, iter) == F(0); ++ipiv)
            ;
        assert(ipiv != md.rows());
        if (ipiv != iter)
            swapRows(md, ipiv, iter);
        divRow(md, iter, md(iter, iter));
        for (size_t i = 0; i < md.rows(); ++i)
            if (i != iter)
                subRows(md, i, iter, md(i, iter));
    }
    
    Matrix<F> minv(m.rows(), m.cols());
    for (size_t i = 0; i < minv.rows(); ++i)
        for (size_t j = 0; j < minv.cols(); ++j)
            minv(i, j) = md(i, j + md.cols() / 2);
    
    assert(minv * m == Matrix<F>::eye(m.rows()));
    return minv;
}

template<class F>
size_t
rank(const Matrix<F>& m)
{
    assert(m.rows() == m.cols());
    Matrix<F> md(m.rows(), 2 * m.cols());
    for (size_t i = 0; i < md.rows(); ++i)
        for (size_t j = 0; j < md.cols() / 2; ++j)
            md(i, j) = m(i, j);
    for (size_t i = 0; i < md.rows(); ++i)
        for (size_t j = md.cols() / 2; j < md.cols(); ++j)
            md(i, j) = 
                (i == j - md.cols() / 2) ? F(1) : F(0);
    for (size_t iter = 0; iter < min(m.rows(), m.cols()); ++iter){
        size_t ipiv = iter;
        for (; ipiv < md.rows() && md(ipiv, iter) == F(0); ++ipiv)
            ;
        if (ipiv == md.rows())
	    return iter;
        if (ipiv != iter)
            swapRows(md, ipiv, iter);
        divRow(md, iter, md(iter, iter));
        for (size_t i = 0; i < md.rows(); ++i)
            if (i != iter)
                subRows(md, i, iter, md(i, iter));
    }
    return m.rows();
}

#endif //MATRIX_H
