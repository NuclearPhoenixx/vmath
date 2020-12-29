/*

*/

#ifndef VMATH_HPP // include guard
#define VMATH_HPP

#include <cmath>

namespace vmath{

    /*
    STANDARD 2D VECTOR FOR VECTOR MATH
    */
    struct Vector2{
        double x,y;

        Vector2(double initx=0., double inity=0.){ //Constructor
            x = initx;
            y = inity;
        }

        Vector2 &operator=(const Vector2 &newvec){ //Assignment Operator
            if(this == &newvec) return *this;
            x = newvec.x;
            y = newvec.y;
            return *this; 
        }

        Vector2 &operator+=(const Vector2 &addvec){ //Simple Vector Addition
            this->x += addvec.x;
            this->y += addvec.y;
            return *this;
        }
        Vector2 &operator-=(const Vector2 &addvec){ //Simple Vector Addition
            this->x -= addvec.x;
            this->y -= addvec.y;
            return *this;
        }
        Vector2 &operator*=(const double num){ //Simple Vector Scalar Multiplication
            this->x *= num;
            this->y *= num;
            return *this;
        }
        Vector2 &operator/=(const double num){ //Simple Vector Scalar Division
            this->x /= num;
            this->y /= num;
            return *this;
        }

        Vector2 operator+(const Vector2 &addvec){ //Simple Vector Addition
            Vector2 newvec = *this;
            newvec += addvec;
            return newvec;
        }
        Vector2 operator-(const Vector2 &addvec){ //Simple Vector Addition
            Vector2 newvec = *this;
            newvec -= addvec;
            return newvec;
        }
        Vector2 operator*(const double num){ //Simple Vector Scalar Multiplication
            Vector2 newvec = *this;
            newvec *= num;
            return newvec;
        }
        Vector2 operator/(const double num){ //Simple Vector Scalar Division
            Vector2 newvec = *this;
            newvec /= num;
            return newvec;
        }

        double length_squared(){ return x*x + y*y; }
        double length(){ return sqrt(length_squared()); }
        
        bool is_normalized(){
            if(this->length_squared() == 1.) return true;
            return false;
        }
        Vector2 normalized(){
            if(this->is_normalized()) return *this;
            Vector2 newvec = *this;
            newvec /= newvec.length();
            return newvec;
        }
        double dot(Vector2 &dotvec){
            double dotprod = 0.;
            dotprod += this->x * dotvec.x;
            dotprod += this->y * dotvec.y;
            return dotprod;
        }
    };

    /*
    STANDARD 3D VECTOR FOR VECTOR MATH
    */
    struct Vector3{
        double x,y,z;

        Vector3(double initx=0., double inity=0.,double initz=0.){ //Constructor
            x = initx;
            y = inity;
            z = initz;
        }

        Vector3 &operator=(const Vector3 &newvec){ //Assignment Operator
            if(this == &newvec) return *this;
            x = newvec.x;
            y = newvec.y;
            z = newvec.z;
            return *this; 
        }

        Vector3 &operator+=(const Vector3 &addvec){ //Simple Vector Addition
            this->x += addvec.x;
            this->y += addvec.y;
            this->z += addvec.z;
            return *this;
        }
        Vector3 &operator-=(const Vector3 &addvec){ //Simple Vector Addition
            this->x -= addvec.x;
            this->y -= addvec.y;
            this->z -= addvec.z;
            return *this;
        }
        Vector3 &operator*=(const double num){ //Simple Vector Scalar Multiplication
            this->x *= num;
            this->y *= num;
            this->z *= num;
            return *this;
        }
        Vector3 &operator/=(const double num){ //Simple Vector Scalar Division
            this->x /= num;
            this->y /= num;
            this->z /= num;
            return *this;
        }

        Vector3 operator+(const Vector3 &addvec){ //Simple Vector Addition
            Vector3 newvec = *this;
            newvec += addvec;
            return newvec;
        }
        Vector3 operator-(const Vector3 &addvec){ //Simple Vector Addition
            Vector3 newvec = *this;
            newvec -= addvec;
            return newvec;
        }
        Vector3 operator*(const double num){ //Simple Vector Scalar Multiplication
            Vector3 newvec = *this;
            newvec *= num;
            return newvec;
        }
        Vector3 operator/(const double num){ //Simple Vector Scalar Division
            Vector3 newvec = *this;
            newvec /= num;
            return newvec;
        }

        double length_squared(){ return x*x + y*y + z*z; }
        double length(){ return sqrt(length_squared()); }
        
        bool is_normalized(){
            if(this->length_squared() == 1.) return true;
            return false;
        }
        Vector3 normalized(){
            if(this->is_normalized()) return *this;
            Vector3 newvec = *this;
            newvec /= newvec.length();
            return newvec;
        }
        double dot(Vector3 &dotvec){
            double dotprod = 0.;
            dotprod += this->x * dotvec.x;
            dotprod += this->y * dotvec.y;
            dotprod += this->z * dotvec.z;
            return dotprod;
        }
        Vector3 cross(Vector3 &crossvec){
            Vector3 &vec = *this;
            double x = vec.y*crossvec.z - vec.z*crossvec.y;
            double y = vec.z*crossvec.x - vec.x*crossvec.z;
            double z = vec.x*crossvec.y - vec.y*crossvec.x;
            return Vector3(x,y,z);
        }
    };

}

#endif
