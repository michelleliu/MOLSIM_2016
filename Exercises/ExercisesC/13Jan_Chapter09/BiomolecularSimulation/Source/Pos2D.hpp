#ifndef _POS_HPP
#define _POS_HPP



#define XDIR 0 
#define YDIR 1 


#define LX 60
#define LY 60


#include <cstring>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>


using namespace std;



class Pos2D{
public:
  Pos2D(int tx, int ty);
  Pos2D();
  Pos2D operator+	(const Pos2D & v) const;
  Pos2D operator-	(const Pos2D & v);
  Pos2D& operator=  (const Pos2D & p);
  bool operator!= (const Pos2D &p)const;
  bool operator== (const Pos2D &p)const;
  void periodicBoundary();
  void periodicSubtraction(Pos2D p1,Pos2D p2);
  void printCout();
  string toString();
  const int operator[] ( const int indx )const;
  int& operator[] ( const int indx );
  // --- data structure ---
  union
  {
    struct
    {
      int x;
      int y;
      
    };
    int	xy[2];
  };
private:
  void copy(const Pos2D& p);
};


// Non-member functions


inline
Pos2D::Pos2D(int tx,int ty){
  x=tx;y=ty;
}

inline
Pos2D::Pos2D(){
  x=0;y=0;
}



inline
Pos2D Pos2D::operator+ (const Pos2D & v)const{
  return Pos2D(x+v.x,y+v.y);
}

inline
Pos2D Pos2D::operator- (const Pos2D & v){
  return Pos2D(x-v.x,y-v.y);
}


inline
bool Pos2D::operator!= (const Pos2D &p)const{
  return (x!=p.x || y!=p.y );
}

inline
bool Pos2D::operator== (const Pos2D &p)const{
  return (x==p.x && y==p.y );
}



inline
const int Pos2D::operator[] (const int indx)const{
  return xy[indx];
}

inline
int &  Pos2D::operator[] (const int indx){
  return xy[indx];
}



inline
Pos2D&	Pos2D::operator= 	( const Pos2D & p )
{
	copy( p );
	return (*this);
}

inline
void Pos2D::copy ( const Pos2D & p ){
	x = p.x;
	y = p.y;
}





inline
void Pos2D::periodicBoundary(){
  if(x>=LX){
    x =x-LX;
  }
  if(x<0){
    x =LX + x;
  }
 if(y>=LY){
    y =y-LY;
  }
  if(y<0){
    y =LY + y;
  }
}

inline
void Pos2D::periodicSubtraction(Pos2D p1,Pos2D p2){
  x=p1.x-p2.x;
  y=p1.y-p2.y;
  if(x > (LX/2)){x=x -LX;}
  else if(x< -(LX/2)){x=x+LX;}
  if(y > (LY/2)){y=y -LY;}
  else if(y< -(LY/2)){y=y+LY;}
}


#endif
