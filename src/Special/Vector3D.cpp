/*
 * =====================================================================================
 *
 *
 *         Author:  Lurker (?), Paul P. Hilscher (2011) 
 *                  from http://www.terathon.com/c4engine/doco/Math/Vector3D.html
 *                  License from post : 
 *                      I'm working on a test 3D engine, and have started on a 3D vector class. 
 *                      I'm looking for opinions, feel free to use it in your code if 
 *                      you wish (probably not enough there to use anyway :p )
 *                  Comment : Nope, very useful ! Thanks :)
 *
 *    Description:  
 *    ToDo       : Add Templates
 * =====================================================================================
 */

#include "Vector3D.h"


Vector3D::Vector3D() {
   x = y = z = 0;
}

Vector3D::Vector3D(double a, double b, double c) {
   x = a;
   y = b;
   z = c;
}

Vector3D::Vector3D(const Vector3D &v) {
   x = v.x;
   y = v.y;
   z = v.z;
}

Vector3D Vector3D::operator+(Vector3D v) {
 return Vector3D((x + v.x), (y + v.y), (z + v.z));
}

Vector3D Vector3D::operator-(Vector3D v) {
 return Vector3D((x - v.x), (y - v.y), (z - v.z));
}

Vector3D Vector3D::operator*(double w) {
 return Vector3D((x * w), (y * w), (z * w));
}

Vector3D Vector3D::operator/(double w) {
 return Vector3D((x / w), (y / w), (z / w));
}

Vector3D Vector3D::operator=(Vector3D v) {
 x = v.x;
 y = v.y;
 z = v.z;
 return *this;
}

void Vector3D::negate() {
 x = -x;
 y = -y;
 z = -z;
}

double Vector3D::length() {
 return sqrt((x * x) + (y * y) + (z * z));
}

void Vector3D::normalize() {
 x /= length();
 y /= length();
 z /= length();
}

void Vector3D::crossWith( Vector3D a )
{
 
    double x1, y1, z1;
 
    x1 = ( a.y * b.z ) - ( a.z * b.y );
    y1 = ( a.z * b.x ) - ( a.x * b.z );
    z1 = ( a.x * b.y ) - ( a.y * b.x );
 
    x = x1;
    y = y1;
    z = z1;
 
}
 
double Vector3D::dotProduct( Vector3D a )
{
 
    return ( x * b.x ) + ( y * b.y ) + ( z * b.z );
 
}

