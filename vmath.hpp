/*
    TODO:
        - More (general) operations for vectors
        - Implement matrices 
*/

#ifndef VMATH_HPP // include guard
#define VMATH_HPP

#include <cmath>
#include <vector>

namespace vmath{

    /*
    STANDARD N-DIMENSIONAL VECTOR FOR VECTOR MATH AND BASE FOR VECTOR2 & VECTOR3
    */
    class VectorN{
        private:
            std::vector<double> vec;

        public:
            VectorN(double n=0.){ //Constructor with dimension
                vec.assign(n,0.);
            }
            VectorN(const std::vector<double> &init){ //Constructor with values
                for(double value:init) vec.push_back(value);
            }

            VectorN &operator=(const VectorN &newvec){ //Assignment Operator
                if(this == &newvec) return *this;
                vec.assign(newvec.dim(),0.);
                for(unsigned int i=0; i<vec.size(); i++) vec[i] = newvec.get(i+1);
                return *this; 
            }

            VectorN &operator+=(const VectorN &addvec){ //Simple Vector Addition
                for(unsigned int i=0; i<vec.size(); i++) vec[i] += addvec.get(i+1);
                return *this;
            }
            VectorN &operator-=(const VectorN &addvec){ //Simple Vector Addition
                for(unsigned int i=0; i<vec.size(); i++) vec[i] -= addvec.get(i+1);
                return *this;
            }
            VectorN &operator*=(const double num){ //Simple Vector Scalar Multiplication
                for(unsigned int i=0; i<vec.size(); i++) vec[i] *= num;
                return *this;
            }
            VectorN &operator/=(const double num){ //Simple Vector Scalar Division
                for(unsigned int i=0; i<vec.size(); i++) vec[i] /= num;
                return *this;
            }

            VectorN operator+(const VectorN &addvec){ //Simple Vector Addition
                VectorN newvec = *this;
                newvec += addvec;
                return newvec;
            }
            VectorN operator-(const VectorN &addvec){ //Simple Vector Addition
                VectorN newvec = *this;
                newvec -= addvec;
                return newvec;
            }
            VectorN operator*(const double num){ //Simple Vector Scalar Multiplication
                VectorN newvec = *this;
                newvec *= num;
                return newvec;
            }
            VectorN operator/(const double num){ //Simple Vector Scalar Division
                VectorN newvec = *this;
                newvec /= num;
                return newvec;
            }

            unsigned int dim() const{ return vec.size(); }
            double get(unsigned int i) const{ return vec[i-1]; }
            void set(unsigned int i, double val) { vec[i] = val; }

            double length_squared(){
                double l = 0.;
                for(double num:vec) l += num*num;
                return l;
            }
            double length(){ return std::sqrt(length_squared()); }
            
            bool is_normalized(){
                if(length_squared() == 1.) return true;
                return false;
            }
            VectorN normalized(){
                if(is_normalized()) return *this;
                VectorN newvec = *this;
                newvec /= newvec.length();
                return newvec;
            }
            double dot(const VectorN &dotvec){
                double dotprod = 0.;
                for(unsigned int i=0; i<vec.size(); i++) dotprod += vec[i] * dotvec.get(i+1);
                return dotprod;
            }
    };

    /*
    STANDARD 2D VECTOR FOR VECTOR MATH
    */
    struct Vector2{
        VectorN vec; //EH, BAD VARIABLE, PLEASE DONT USE LOL

        Vector2(const double initx=0., const double inity=0.){ //Constructor
            //std::vector<double> v = {initx,inity};
            vec = VectorN({initx,inity});
        }
        Vector2(const VectorN v){ //"Copy" Constructor
            vec = v;
        }

        Vector2 &operator=(const Vector2 &newvec){ //Assignment Operator
            if(this == &newvec) return *this;
            vec = newvec.vec;
            return *this; 
        }
        Vector2 &operator+=(const Vector2 &addvec){ //Simple Vector Addition
            vec += addvec.vec;
            return *this;
        }
        Vector2 &operator-=(const Vector2 &addvec){ //Simple Vector Addition
            vec -= addvec.vec;
            return *this;
        }
        Vector2 &operator*=(const double num){ //Simple Vector Scalar Multiplication
            vec *= num;
            return *this;
        }
        Vector2 &operator/=(const double num){ //Simple Vector Scalar Division
            vec /= num;
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

        double x(){ return vec.get(1); }
        double y(){ return vec.get(2); }

        double length_squared(){ return vec.length_squared(); }
        double length(){ return vec.length(); }

        bool is_normalized(){ return vec.is_normalized(); }
        Vector2 normalized(){ return Vector2(vec.normalized()); }
        double dot(const Vector2 &dotvec){ return vec.dot(dotvec.vec); }
    };

    /*
    STANDARD 3D VECTOR FOR VECTOR MATH
    */
    struct Vector3{
        VectorN vec; //EH, BAD VARIABLE, PLEASE DONT USE LOL

        Vector3(const double initx=0., const double inity=0., const double initz=0.){ //Constructor
            //std::vector<double> v = {initx,inity,initz};
            vec = VectorN({initx,inity,initz});
        }
        Vector3(const VectorN v){ //"Copy" Constructor
            vec = v;
        }

        Vector3 &operator=(const Vector3 &newvec){ //Assignment Operator
            if(this == &newvec) return *this;
            vec = newvec.vec;
            return *this; 
        }
        Vector3 &operator+=(const Vector3 &addvec){ //Simple Vector Addition
            vec += addvec.vec;
            return *this;
        }
        Vector3 &operator-=(const Vector3 &addvec){ //Simple Vector Addition
            vec -= addvec.vec;
            return *this;
        }
        Vector3 &operator*=(const double num){ //Simple Vector Scalar Multiplication
            vec *= num;
            return *this;
        }
        Vector3 &operator/=(const double num){ //Simple Vector Scalar Division
            vec /= num;
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

        double x(){ return vec.get(1); }
        double y(){ return vec.get(2); }
        double z(){ return vec.get(3); }

        double length_squared(){ return vec.length_squared(); }
        double length(){ return vec.length(); }

        bool is_normalized(){ return vec.is_normalized(); }
        Vector3 normalized(){ return Vector3(vec.normalized()); }
        double dot(const Vector3 &dotvec){ return vec.dot(dotvec.vec); }

        Vector3 cross(Vector3 &crossvec){
            double newx = y()*crossvec.z() - z()*crossvec.y();
            double newy = z()*crossvec.x() - x()*crossvec.z();
            double newz = x()*crossvec.y() - y()*crossvec.x();
            return Vector3(newx,newy,newz);
        }
    };

}

#endif
