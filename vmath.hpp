/*
    TODO:
        - More (general) operations for vectors and matrices
        - Inheritance and (some) templates would be MUCH cleaner
        - Matrix2 and Matrix3 with optional inverse functions
*/

#ifndef VMATH_HPP // include guard
#define VMATH_HPP

#include <cmath> //used only TWICE, once for std::pow and once std::sqrt....
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
                vec = init;
            }

            VectorN &operator=(const VectorN &newvec){ //Assignment Operator
                if(this == &newvec) return *this;
                vec = newvec.vec;
                //vec.assign(newvec.dim(),0.);
                //for(size_t i=0; i<vec.size(); i++) vec[i] = newvec.get(i+1);
                return *this; 
            }

            VectorN &operator+=(const VectorN &addvec){ //Simple Vector Addition
                for(size_t i=0; i<vec.size(); i++) vec[i] += addvec.get(i+1);
                return *this;
            }
            VectorN &operator-=(const VectorN &addvec){ //Simple Vector Addition
                for(size_t i=0; i<vec.size(); i++) vec[i] -= addvec.get(i+1);
                return *this;
            }
            VectorN &operator*=(const double num){ //Simple Vector Scalar Multiplication
                for(size_t i=0; i<vec.size(); i++) vec[i] *= num;
                return *this;
            }
            VectorN &operator/=(const double num){ //Simple Vector Scalar Division
                for(size_t i=0; i<vec.size(); i++) vec[i] /= num;
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
            void set(unsigned int i, double val) { vec[i-1] = val; }

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
                for(size_t i=0; i<vec.size(); i++) dotprod += vec[i] * dotvec.get(i+1);
                return dotprod;
            }
    };

    /*
    STANDARD 2D VECTOR FOR VECTOR MATH
    */
    struct Vector2{
        VectorN vec; //EH, BAD VARIABLE, PLEASE DONT USE LOL

        Vector2(const double initx=0., const double inity=0.){ //Constructor
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

        void set(unsigned int i, double val) { vec.set(i,val); }
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

    class MatrixN{
        private:
            std::vector<std::vector<double>> mat;
            
            void remove_line(unsigned int i=0){ mat.erase(mat.begin()+i-1); }
            void remove_column(unsigned int i=0){ for(auto &matrix:mat) matrix.erase(matrix.begin()+i-1); }

            double subdet(unsigned int i, unsigned int j, MatrixN matrix){
                matrix.remove_line(i);
                matrix.remove_column(j);
                return matrix.det();
            }

        public:
            MatrixN(double n=0.,double m=0.){ //Constructor with dimension
                mat.assign(n,std::vector<double>(m,0.));
            }
            MatrixN(const std::vector<std::vector<double>> &init){ //Constructor with values
                mat = init;
            }
            
            MatrixN &operator=(const MatrixN &newmat){
                if(this == &newmat) return *this;
                mat.assign((newmat.dim()).x(),std::vector<double>((newmat.dim()).y(),0.));
                for(size_t n=0; n<mat.size(); n++){
                    for(size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] = newmat.get(n+1,m+1);
                    }
                }
                return *this; 
            }
            
            MatrixN &operator+=(const MatrixN &addmat){
                for(size_t n=0; n<mat.size(); n++){
                    for(size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] += addmat.get(n+1, m+1);
                    }
                }
                return *this;
            }
            MatrixN &operator-=(const MatrixN &addmat){
                for(size_t n=0; n<mat.size(); n++){
                    for(size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] -= addmat.get(n+1, m+1);
                    }
                }
                return *this;
            }
            MatrixN &operator*=(const double num){
                for(size_t n=0; n<mat.size(); n++){
                    for(size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] *= num;
                    }
                }
                return *this;
            }
            MatrixN &operator*=(const MatrixN &multmat){
                MatrixN tempmat = mat;
                for(size_t i=0; i<mat.size(); i++){
                    for(size_t k=0; k<mat[i].size(); k++){
                        double val = 0.;
                        for(size_t j=0; j<mat[i].size(); j++){
                            val += mat[i][j] * multmat.get(j+1,k+1);
                        }
                        tempmat.set(i+1,k+1,val);
                    }
                }
                *this = tempmat; //meh a bit ugly
                return *this;
            }
            MatrixN &operator/=(const double num){
                for(size_t n=0; n<mat.size(); n++){
                    for(size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] /= num;
                    }
                }
                return *this;
            }
            
            MatrixN operator+(const MatrixN &addmat){
                MatrixN newvec = *this;
                newvec += addmat;
                return newvec;
            }
            MatrixN operator-(const MatrixN &addmat){
                MatrixN newvec = *this;
                newvec -= addmat;
                return newvec;
            }
            MatrixN operator*(const double num){ //Matrix Scalar Operation
                MatrixN newvec = *this;
                newvec *= num;
                return newvec;
            }
            MatrixN operator*(const MatrixN &multmat){ //Matrix Matrix Operation
                MatrixN newvec = *this;
                newvec *= multmat;
                return newvec;
            }
            VectorN operator*(const VectorN &multvec){ //Matrix Vector Operation
                VectorN newvec(multvec.dim());
                for(size_t i=1; i<=multvec.dim(); i++){
                    double val = 0.;
                    for(size_t j=1; j<=dim().y(); j++){
                        val += get(i,j) * multvec.get(j);
                    }
                    newvec.set(i,val);
                }
                return newvec;
            }
            MatrixN operator/(const double num){
                MatrixN newvec = *this;
                newvec /= num;
                return newvec;
            }

            Vector2 dim() const{ return Vector2(mat.size(),mat[0].size()); }
            double get(unsigned int n, unsigned int m) const{ return mat[n-1][m-1]; }
            void set(unsigned int n, unsigned int m, double val) { mat[n-1][m-1] = val; }

            double det(unsigned int i=1, unsigned int j=1){
                if(dim().x() == 2) return get(1,1) * get(2,2) - get(2,1) * get(1,2);

                double val = 0.;
                for(size_t jj=j; jj<=(dim()).y(); jj++){
                    val += std::pow(-1,i+jj) * get(i,jj) * subdet(i,jj,*this);
                }
                return val;
            }
    };

}

#endif
