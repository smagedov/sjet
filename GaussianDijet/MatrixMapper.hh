#ifndef MATRIXMAPPER_HH_
#define MATRIXMAPPER_HH_

#include <cassert>
#include <limits>

class MatrixMapper
{
public:
    inline MatrixMapper(const unsigned i_nRows, const unsigned i_nCols)
        : nRows_(i_nRows), nCols_(i_nCols)
    {
        assert(nRows_);
        assert(nCols_);
        const unsigned long sz = static_cast<unsigned long>(nRows_)*nCols_;
        assert(sz <= std::numeric_limits<unsigned>::max());
        size_ = sz;
    }

    inline unsigned nRows() const {return nRows_;}
    inline unsigned nCols() const {return nCols_;}
    inline unsigned size() const {return size_;}

    inline unsigned operator()(const unsigned row, const unsigned col) const
    {
        assert(row < nRows_);
        assert(col < nCols_);
        return row*nCols_ + col;
    }

    inline void inverse(const unsigned linearIndex,
                        unsigned* row, unsigned* col) const
    {
        assert(linearIndex < size_);
        if (row)
            *row = linearIndex / nCols_;
        if (col)
            *col = linearIndex % nCols_;
    }

private:
    unsigned nRows_;
    unsigned nCols_;
    unsigned size_;
};

#endif // MATRIXMAPPER_HH_
