#pragma once
#include <vector>
#include <iostream>

template< typename ValueType = double,
	  typename IndexType = int,
	  typename ValueContainer = std::vector< ValueType >,
	  typename IndexContainer = std::vector< IndexType > >
class CSR
{
public:
   CSR () {}
   CSR ( int _nRow, int _nCol, int _nNZV )
      : mRow( _nRow ), mCol( _nCol ), mNNZ( _nNZV )
   {}
   
   IndexContainer& Row ()  { return mVecRow; }
   IndexContainer& Col ()  { return mVecCol; }
   ValueContainer& NZV ()  { return mVecNZV;  }

   const IndexContainer& Row () const  { return mVecRow; }
   const IndexContainer& Col () const  { return mVecCol; }
   const ValueContainer& NZV () const  { return mVecNZV; }

   int& nRow ()    { return mRow; }
   int& nCol ()    { return mCol; }
   int& nNZV ()    { return mNNZ; }
   
   const int nRow () const   { return mRow; }
   const int nCol () const   { return mCol; }
   const int nNZV () const   { return mNNZ; }
   
   friend std::ostream& operator << ( std::ostream& os, const CSR<>& csr );

private:
   int mRow;
   int mCol;
   int mNNZ;

   IndexContainer mVecRow;
   IndexContainer mVecCol;
   ValueContainer mVecNZV;
};

std::ostream& operator<< ( std::ostream& os, const CSR<>& csr )	  
{
   os << "+++++++++++++++ CSR Information +++++++++++++++++" << std::endl;
   os << "nRow = " << csr.nRow() << std::endl;
   os << "nCol = " << csr.nCol() << std::endl;
   os << "nNZV = " << csr.nNZV() << std::endl;

   os << "Row Content --- --- --- " << std::endl;
   for( std::size_t i = 0; i < csr.Row().size(); ++i )
   {
      os << i << "\t" << csr.Row()[i] << std::endl;
   }
   os << "Column Content --- --- --- " << std::endl;
   for( std::size_t i = 0; i < csr.Col().size(); ++i )
   {
      os << i << "\t" << csr.Col()[i] << std::endl;
   }

   os << "NZV Content --- --- --- " << std::endl;
   for( std::size_t i = 0; i < csr.NZV().size(); ++i )
   {
      os << i << "\t" << csr.NZV()[i] << std::endl;
   }

   os << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;

   return os;
}
   
