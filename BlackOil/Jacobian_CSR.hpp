#pragma once
#include <vector>

template< typename ValueType = double,
	  typename IndexType = int,
	  typename ValueContainer = std::vector< ValueType >,
	  typename IndexContainer = std::vector< IndexType > >
class CSR
{
public:
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
   
   // void Update_Status ()
   // {
   //    mRow = mVecRow.size();
   //    mCol = mVecCol.size();
   //    mNNZ = mVecNZV.size();   
   // }
   
private:
   int mRow;
   int mCol;
   int mNNZ;

   IndexContainer mVecRow;
   IndexContainer mVecCol;
   ValueContainer mVecNZV;
};
	  
