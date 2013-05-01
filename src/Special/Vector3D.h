/*
* =====================================================================================
*
*
*         Author:  Lurker (?), Paul P. Hilscher (2011)
*
*                  Using from
*                  from http://www.terathon.com/c4engine/doco/Math/Vector3D.html
*                  License from post : 
*                      I'm working on a test 3D engine, and have started on a 3D vector class. 
*                      I'm looking for opinions, feel free to use it in your code if 
*                      you wish (probably not enough there to use anyway :p )
*
*    Description:  
*    ToDo       : Add Templates
* =====================================================================================
*/


#ifndef __VECTOR3D_H
#define __VECTOR3D_H

/**
*   @brief a three dimensional vector class for Vector3D class
*
**/
class Vector3D {

 public:

  double x, ///< x-component
         y, ///< y-component
         z; ///< z-component

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D() {
   x = y = z = 0;
  }

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D(double a, double b, double c) 
  {

   x = a;
   y = b;
   z = c;

  }

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D(const Vector3D &v) 
  {

   x = v.x;
   y = v.y;
   z = v.z;

  }

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D operator+(Vector3D v);

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D operator-(Vector3D v);

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D operator*(double w);

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D operator/(double w);

  /**
  *
  * @brief Constructors
  *
  **/
  Vector3D operator=(Vector3D v);

  /**
  *
  * @brief Constructors
  *
  **/
  void negate();
  
  /**
  *
  * @brief Constructors
  *
  **/
  double length();
  
  /**
  * @brief Constructors
  **/
  void normalize();

  /**
  * @brief Constructors
  **/
  void crossWith( Vector3D a );

  /**
  * @brief Constructors
  **/
  double dotProduct( Vector3D a );
};


#endif // VECTOR_3D
