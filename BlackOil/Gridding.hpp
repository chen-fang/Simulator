#pragma once
#include <cstddef>

#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADscalar;
typedef adetl::ADvector   ADvector;

struct Cartesian3D
{
   Cartesian3D ( std::size_t _dx, std::size_t _dy, std::size_t _dz,
		 std::size_t _nx, std::size_t _ny, std::size_t _nz )
      : Dx(_dx), Dy(_dy), Dz(_dz), Nx(_nx), Ny(_ny), Nz(_nz)
   {
      std::cout << "Cartesian 3D constructor" << std::endl;
   }
   
   void Set( std::size_t _dx, std::size_t _dy, std::size_t _dz,
	     std::size_t _nx, std::size_t _ny, std::size_t _nz )
   {
      Dx = _dx;
      Dy = _dy;
      Dz = _dz;
      Nx = _nx;
      Ny = _ny;
      Nz = _nz;
   }

   std::size_t GetGridBlockIndex ( std::size_t _i, std::size_t _j, std::size_t _k ) const
   {
      return _k*Ny*Nx + _j*Nx + _i;
   }

   std::size_t GetGridBlockNumber () const
   {
      return Nx * Ny * Nz;
   }

   std::size_t Dx;
   std::size_t Dy;
   std::size_t Dz;

   std::size_t Nx;
   std::size_t Ny;
   std::size_t Nz;
};


struct ConnectionList
{
   ConnectionList ( const Cartesian3D& _grid )
   {
      Initialize( _grid );
   }
   
   void Initialize ( const Cartesian3D& _grid );
   
   std::size_t Size;
   std::size_t XSize;
   std::size_t YSize;
   std::size_t ZSize;
   
   struct pair
   {
      std::size_t left;
      std::size_t right;
   };
   std::vector< pair > CList;

   const pair& operator[] ( std::size_t i )       { return CList[i]; }
   const pair& operator[] ( std::size_t i ) const { return CList[i]; }

};

void ConnectionList :: Initialize( const Cartesian3D& _grid )
{
   std::size_t _Nx = _grid.Nx;
   std::size_t _Ny = _grid.Ny;
   std::size_t _Nz = _grid.Nz;

   Size  = 3*_Nx*_Ny*_Nz - _Nx*_Ny - _Nx*_Nz - _Ny*_Nz;
   XSize = (_Nx-1)*_Ny*_Nz;;
   YSize = (_Ny-1)*_Nx*_Nz;;
   ZSize = (_Nz-1)*_Nx*_Ny;;

   CList.resize ( Size );
   std::size_t index = 0;

   // 1st direction in all layers
   for (std::size_t k=0; k < _Nz; ++k)
   {
      for (std::size_t j=0; j < _Ny; ++j)
      {
	 for (std::size_t i=0; i < _Nx-1; ++i)
	 {
	    CList[index].left = i + j*_Nx + k*_Nx*_Ny;
	    CList[index].right = CList[index].left + 1;
	    index++;
	 }
      }
   }
   // 2nd direction in all layers
   for (std::size_t k=0; k < _Nz; ++k)
   {
      for (std::size_t j=0; j < _Ny-1; ++j)
      {
	 for (std::size_t i=0; i < _Nx; ++i)
	 {
	    CList[index].left = i + j*_Nx + k*_Nx*_Ny;
	    CList[index].right = CList[index].left + _Nx;
	    index++;
	 }
      }
   }
   // 3rd direction
   for (std::size_t k=0; k < _Nz-1; ++k)
   {
      for (std::size_t j=0; j < _Ny; ++j)
      {
	 for (std::size_t i=0; i < _Nx; ++i)
	 {
	    CList[index].left = i + j*_Nx + k*_Nx*_Ny;
	    CList[index].right = CList[index].left + _Nx*_Ny;
	    index++;
	 }
      }
   }
}
