/*
    vmath.hpp

    Simple to use vector and matrix math library.

    2021, Phoenix1747, MIT License.
    
    Github: https://github.com/Phoenix1747/vmath

    (TODO) Inheritance and (some) templates would be MUCH cleaner
        - Orthonormalization -> Gram-Schmidt
*/

#ifndef VMATH_HPP
#define VMATH_HPP

#include <cmath> //used only for pow, abs, sqrt, sin, cos, acos, M_PI
#include <vector>

namespace vmath{

    /*
    STANDARD N-DIMENSIONAL VECTOR FOR VECTOR MATH AND BASE FOR VECTOR2 & VECTOR3
    */
    class VectorN{
        //private:
            //std::vector<double> vec;

        public:
            std::vector<double> vec;

            VectorN(double n=0.){ //Constructor for vector dimension n
                vec.assign(n,0.);
            }
            VectorN(const std::vector<double> &init){
                vec = init;
            }

            VectorN &operator=(const VectorN &newvec){
                if(this == &newvec) return *this;
                vec = newvec.vec;
                //vec.assign(newvec.dim(),0.);
                //for(std::size_t i=0; i<vec.size(); i++) vec[i] = newvec.get(i+1);
                return *this; 
            }
            bool operator==(const VectorN &vector){
                if(vec != vector.vec) return false;
                return true;
            }
            bool operator!=(const VectorN &vector){
                if(vec == vector.vec) return false;
                return true;
            }

            VectorN &operator+=(const VectorN &addvec){ //Vector Addition
                for(std::size_t i=0; i<vec.size(); i++) vec[i] += addvec.get(i+1);
                return *this;
            }
            VectorN &operator-=(const VectorN &addvec){ //Vector Subtraction
                for(std::size_t i=0; i<vec.size(); i++) vec[i] -= addvec.get(i+1);
                return *this;
            }
            VectorN &operator*=(const double num){ //Vector Scalar Multiplication
                for(std::size_t i=0; i<vec.size(); i++) vec[i] *= num;
                return *this;
            }
            VectorN &operator/=(const double num){ //Vector Scalar Division
                for(std::size_t i=0; i<vec.size(); i++) vec[i] /= num;
                return *this;
            }

            VectorN operator+(const VectorN &addvec){
                VectorN newvec = *this;
                newvec += addvec;
                return newvec;
            }
            VectorN operator-(const VectorN &addvec){
                VectorN newvec = *this;
                newvec -= addvec;
                return newvec;
            }
            VectorN operator*(const double num){
                VectorN newvec = *this;
                newvec *= num;
                return newvec;
            }
            VectorN operator/(const double num){
                VectorN newvec = *this;
                newvec /= num;
                return newvec;
            }

            unsigned int dim() const{ return vec.size(); } //return vector dimension
            double get(unsigned int i) const{ return vec[i-1]; }
            void set(unsigned int i, double val) { vec[i-1] = val; }

            double length_squared() const{
                double l = 0.;
                for(double num:vec) l += num*num;
                return l;
            }
            double length() const{ return std::sqrt(length_squared()); }
            
            bool is_normalized(){ //if length() == 1
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
                for(std::size_t i=0; i<vec.size(); i++) dotprod += vec[i] * dotvec.get(i+1);
                return dotprod;
            }
            VectorN abs(){
                std::vector<double> newvec;
                for(auto val:vec) newvec.push_back(std::abs(val));
                return VectorN(newvec);
            }
            double distance_squared_to(const VectorN &tovec){ return (*this-tovec).length_squared(); }
            double distance_to(const VectorN &tovec){ return (*this-tovec).length(); }
    };

    /*
    STANDARD 2D VECTOR FOR VECTOR MATH
    */
    struct Vector2{
        VectorN vec;

        Vector2(const double initx=0., const double inity=0.){ //Create vector with two coords
            vec = VectorN({initx,inity});
        }
        Vector2(const VectorN &v){ //Create vector out of existing VectorN
            vec = v;
        }

        Vector2 &operator=(const Vector2 &newvec){
            if(this == &newvec) return *this;
            vec = newvec.vec;
            return *this; 
        }
        bool operator==(const Vector2 &vector){
            //if(vec != vector.vec) return false;
            //return true;
            return (vec == vector.vec);
        }
        bool operator!=(const Vector2 &vector){
            //if(vec == vector.vec) return false;
            //return true;
            return (vec != vector.vec);
        }

        Vector2 &operator+=(const Vector2 &addvec){
            vec += addvec.vec;
            return *this;
        }
        Vector2 &operator-=(const Vector2 &addvec){
            vec -= addvec.vec;
            return *this;
        }
        Vector2 &operator*=(const double num){
            vec *= num;
            return *this;
        }
        Vector2 &operator/=(const double num){
            vec /= num;
            return *this;
        }

        Vector2 operator+(const Vector2 &addvec){
            Vector2 newvec = *this;
            newvec += addvec;
            return newvec;
        }
        Vector2 operator-(const Vector2 &addvec){
            Vector2 newvec = *this;
            newvec -= addvec;
            return newvec;
        }
        Vector2 operator*(const double num){
            Vector2 newvec = *this;
            newvec *= num;
            return newvec;
        }
        Vector2 operator/(const double num){
            Vector2 newvec = *this;
            newvec /= num;
            return newvec;
        }

        double x(){ return vec.get(1); } //return first coord
        double y(){ return vec.get(2); } //return second coord

        double length_squared(){ return vec.length_squared(); }
        double length() const{ return vec.length(); }

        void set(unsigned int i, double val) { vec.set(i,val); }
        bool is_normalized(){ return vec.is_normalized(); }
        Vector2 normalized(){ return Vector2(vec.normalized()); }
        double dot(const Vector2 &dotvec){ return vec.dot(dotvec.vec); }
        Vector2 abs(){ return Vector2(vec.abs()); }

        double angle_to(const Vector2 &tovec){
            const double dotprod = dot(tovec);
            const double lengthprod = length() * tovec.length();
            return std::acos(dotprod/lengthprod);
        }
        double aspect(){ return x()/y(); }
        Vector2 rotated(double phi){
            const double xi = x() * std::cos(phi) - y() * std::sin(phi);
            const double yi = y() * std::cos(phi) + x() * std::sin(phi);
            return Vector2(xi,yi);
        }
        Vector2 tangent(){ return rotated(M_PI/2.); }
        double distance_squared_to(const Vector2 &tovec){ return vec.distance_squared_to(tovec.vec); }
        double distance_to(const Vector2 &tovec){ return vec.distance_to(tovec.vec); }
        
    };

    /*
    STANDARD 3D VECTOR FOR VECTOR MATH
    */
    struct Vector3{
        VectorN vec;

        Vector3(const double initx=0., const double inity=0., const double initz=0.){ //Create vector with the 3 coords
            vec = VectorN({initx,inity,initz});
        }
        Vector3(const VectorN &v){ //Create vector out of existing vector
            vec = v;
        }

        Vector3 &operator=(const Vector3 &newvec){
            if(this == &newvec) return *this;
            vec = newvec.vec;
            return *this;
        }
        bool operator==(const Vector3 &vector){
            //if(vec != vector.vec) return false;
            //return true;
            return (vec == vector.vec);
        }
        bool operator!=(const Vector3 &vector){
            //if(vec == vector.vec) return false;
            //return true;
            return (vec != vector.vec);
        }

        Vector3 &operator+=(const Vector3 &addvec){
            vec += addvec.vec;
            return *this;
        }
        Vector3 &operator-=(const Vector3 &addvec){
            vec -= addvec.vec;
            return *this;
        }
        Vector3 &operator*=(const double num){
            vec *= num;
            return *this;
        }
        Vector3 &operator/=(const double num){
            vec /= num;
            return *this;
        }

        Vector3 operator+(const Vector3 &addvec){
            Vector3 newvec = *this;
            newvec += addvec;
            return newvec;
        }
        Vector3 operator-(const Vector3 &addvec){
            Vector3 newvec = *this;
            newvec -= addvec;
            return newvec;
        }
        Vector3 operator*(const double num){
            Vector3 newvec = *this;
            newvec *= num;
            return newvec;
        }
        Vector3 operator/(const double num){
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
        Vector3 abs(){ return Vector3(vec.abs()); }

        Vector3 cross(Vector3 &crossvec){
            const double newx = y()*crossvec.z() - z()*crossvec.y();
            const double newy = z()*crossvec.x() - x()*crossvec.z();
            const double newz = x()*crossvec.y() - y()*crossvec.x();
            return Vector3(newx,newy,newz);
        }
        double distance_squared_to(const Vector3 &tovec){ return vec.distance_squared_to(tovec.vec); }
        double distance_to(const Vector3 &tovec){ return vec.distance_to(tovec.vec); }

    };

    /*
    STANDARD N-DIMENSIONAL MATRIX FOR VECTOR/MATRIX MATH AND BASE FOR MATRIX2 & MATRIX3
    */
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
            MatrixN(double n=0.,double m=0.){ //Create (n x m) matrix with values 0
                mat.assign(n,std::vector<double>(m,0.));
            }
            MatrixN(const std::vector<std::vector<double>> &init){ //Create matrix out of existing 2-dimensional vector
                mat = init;
            }
            
            MatrixN &operator=(const MatrixN &newmat){
                if(this == &newmat) return *this;
                mat.assign((newmat.dim()).x(),std::vector<double>((newmat.dim()).y(),0.));
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] = newmat.get(n+1,m+1);
                    }
                }
                return *this; 
            }
            bool operator==(const MatrixN &matrix){
                if(mat != matrix.mat) return false;
                return true;
            }
            bool operator!=(const MatrixN &matrix){
                if(mat == matrix.mat) return false;
                return true;
            }
            
            MatrixN &operator+=(const MatrixN &addmat){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] += addmat.get(n+1, m+1);
                    }
                }
                return *this;
            }
            MatrixN &operator-=(const MatrixN &addmat){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] -= addmat.get(n+1, m+1);
                    }
                }
                return *this;
            }
            MatrixN &operator*=(const double num){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat[n].size(); m++){
                        mat[n][m] *= num;
                    }
                }
                return *this;
            }
            MatrixN &operator*=(const MatrixN &multmat){
                MatrixN tempmat = mat;
                for(std::size_t i=0; i<mat.size(); i++){
                    for(std::size_t k=0; k<mat[i].size(); k++){
                        double val = 0.;
                        for(std::size_t j=0; j<mat[i].size(); j++){
                            val += mat[i][j] * multmat.get(j+1,k+1);
                        }
                        tempmat.set(i+1,k+1,val);
                    }
                }
                *this = tempmat; //meh a bit ugly
                return *this;
            }
            MatrixN &operator/=(const double num){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat[n].size(); m++){
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
                for(std::size_t i=1; i<=multvec.dim(); i++){
                    double val = 0.;
                    for(std::size_t j=1; j<=dim().y(); j++){
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
                for(std::size_t jj=j; jj<=(dim()).y(); jj++){
                    val += std::pow(-1,i+jj) * get(i,jj) * subdet(i,jj,*this);
                }
                return val;
            }
            MatrixN transpose(){
                MatrixN newmat = *this;
                for(std::size_t n=1; n<=dim().x(); n++){
                    for(std::size_t m=1; m<=dim().y(); m++) newmat.set(n,m,get(m,n));
                }
                return newmat;
            }
            bool is_normalized(){ //if length of every column == 1
                for(std::size_t m=1; m<=dim().y(); m++){
                    std::vector<double> tempvec;
                    for(std::size_t n=1; n<=dim().x(); n++) tempvec.push_back(get(n,m));
                    VectorN vector(tempvec);
                    if(!vector.is_normalized()) return false;
                }
                return true;
            }
            bool is_orthogonalized(){ //if every column is perpendicular to every other one
                std::vector<VectorN> vecholder;
                
                for(std::size_t m=1; m<=dim().y(); m++){
                    std::vector<double> tempvec;
                    for(std::size_t n=1; n<=dim().x(); n++) tempvec.push_back(get(n,m));
                    vecholder.push_back(tempvec);
                }
                for(auto vector:vecholder){
                    for(auto anothervec:vecholder){
                        if(vector == anothervec) continue;
                        if(vector.dot(anothervec) != 0) return false;
                    }
                }
                return true;
            }
            bool is_orthonormalized(){
                //if(is_orthogonalized() && is_normalized()) return true;
                //return false;
                return (is_orthogonalized() && is_normalized());
            }
    };

    /*
    STANDARD 2x2 MATRIX FOR VECTOR MATH
    */
    struct Matrix2{
        MatrixN mat;

        Matrix2(double a11=0., double a12=0., double a21=0., double a22=0.){ //Construct with matrix values
            mat = MatrixN({{a11,a12},{a21,a22}});
        }
        Matrix2(const std::vector<std::vector<double>> &init){ //Construct with 2d 2x2 vector
            mat = MatrixN(init);
        }
        Matrix2(const MatrixN &init){ //Construct with existing 2x2 MatrixN
            mat = init;
        }
            
        Matrix2 &operator=(const Matrix2 &newmat){
            if(this == &newmat) return *this;
            mat = newmat.mat;
            return *this; 
        }
        bool operator==(const Matrix2 &matrix){
            return (mat == matrix.mat);
        }
        bool operator!=(const Matrix2 &matrix){
            return (mat != matrix.mat);
        }
            
        Matrix2 &operator+=(const Matrix2 &addmat){
            mat += addmat.mat;
            return *this;
        }
        Matrix2 &operator-=(const Matrix2 &addmat){
            mat -= addmat.mat;
            return *this;
        }
        Matrix2 &operator*=(const double num){
            mat *= num;
            return *this;
        }
        Matrix2 &operator*=(const Matrix2 &multmat){
            mat *= multmat.mat;
            return *this;
        }
        Matrix2 &operator/=(const double num){
            mat /= num;
            return *this;
        }
        
        Matrix2 operator+(const Matrix2 &addmat){
            Matrix2 newmat = *this;
            newmat += addmat;
            return newmat;
        }
        Matrix2 operator-(const Matrix2 &addmat){
            Matrix2 newmat = *this;
            newmat -= addmat;
            return newmat;
        }
        Matrix2 operator*(const double num){ //Matrix Scalar Operation
            Matrix2 newmat = *this;
            newmat *= num;
            return newmat;
        }
        Matrix2 operator*(const Matrix2 &multmat){ //Matrix Matrix Operation
            Matrix2 newvec = *this;
            newvec *= multmat;
            return newvec;
        }
        Vector2 operator*(const Vector2 &multvec){ //Matrix Vector Operation
            return Vector2(mat * multvec.vec);;
        }
        Matrix2 operator/(const double num){
            Matrix2 newmat = *this;
            newmat /= num;
            return newmat;
        }

        Vector2 dim() const{ return mat.dim(); }
        double get(unsigned int n, unsigned int m) const{ return mat.get(n,m); }
        void set(unsigned int n, unsigned int m, double val) { mat.set(n,m,val); }
        double det(){ return mat.det(); }
        Matrix2 transpose(){ return Matrix2(mat.transpose()); }
        bool is_normalized(){ return mat.is_normalized(); }
        bool is_orthogonalized(){ return mat.is_orthogonalized(); }
        bool is_orthonormalized(){ return mat.is_orthonormalized(); }
        Matrix2 inverse(){
            Matrix2 newmat;
            newmat.set(1,1,get(2,2));
            newmat.set(1,2,-get(1,2));
            newmat.set(2,1,-get(2,1));
            newmat.set(2,2,get(1,1));
            newmat /= det();
            return newmat;            
        }

    };

    /*
    STANDARD 3x3 MATRIX FOR VECTOR MATH
    */
    struct Matrix3{
        MatrixN mat;

        //Construct with matrix values
        Matrix3(double a11=0.,double a12=0.,double a13=0.,double a21=0.,double a22=0.,double a23=0.,double a31=0.,double a32=0.,double a33=0.){
            mat = MatrixN({{a11,a12,a13},{a21,a22,a23},{a31,a32,a33}});
        }
        Matrix3(const std::vector<std::vector<double>> &init){ //Construct with 2d 3x3 vector
            mat = MatrixN(init);
        }
        Matrix3(const MatrixN &init){ //Construct with existing 3x3 MatrixN
            mat = init;
        }
            
        Matrix3 &operator=(const Matrix3 &newmat){
            if(this == &newmat) return *this;
            mat = newmat.mat;
            return *this;
        }
        bool operator==(const Matrix3 &matrix){
            return (mat == matrix.mat);
        }
        bool operator!=(const Matrix3 &matrix){
            return (mat != matrix.mat);
        }
        
        Matrix3 &operator+=(const Matrix3 &addmat){
            mat += addmat.mat;
            return *this;
        }
        Matrix3 &operator-=(const Matrix3 &addmat){
            mat -= addmat.mat;
            return *this;
        }
        Matrix3 &operator*=(const double num){
            mat *= num;
            return *this;
        }
        Matrix3 &operator*=(const Matrix3 &multmat){
            mat *= multmat.mat;
            return *this;
        }
        Matrix3 &operator/=(const double num){
            mat /= num;
            return *this;
        }
        
        Matrix3 operator+(const Matrix3 &addmat){
            Matrix3 newmat = *this;
            newmat += addmat;
            return newmat;
        }
        Matrix3 operator-(const Matrix3 &addmat){
            Matrix3 newmat = *this;
            newmat -= addmat;
            return newmat;
        }
        Matrix3 operator*(const double num){ //Matrix Scalar Operation
            Matrix3 newmat = *this;
            newmat *= num;
            return newmat;
        }
        Matrix3 operator*(const Matrix3 &multmat){ //Matrix Matrix Operation
            Matrix3 newvec = *this;
            newvec *= multmat;
            return newvec;
        }
        Vector3 operator*(const Vector3 &multvec){ //Matrix Vector Operation
            return Vector3(mat * multvec.vec);;
        }
        Matrix3 operator/(const double num){
            Matrix3 newmat = *this;
            newmat /= num;
            return newmat;
        }

        Vector2 dim() const{ return mat.dim(); }
        double get(unsigned int n, unsigned int m) const{ return mat.get(n,m); }
        void set(unsigned int n, unsigned int m, double val) { mat.set(n,m,val); }
        double det(){ return mat.det(); }
        Matrix2 transpose(){ return Matrix2(mat.transpose()); }
        bool is_normalized(){ return mat.is_normalized(); }
        bool is_orthogonalized(){ return mat.is_orthogonalized(); }
        bool is_orthonormalized(){ return mat.is_orthonormalized(); }
        Matrix3 inverse(){
            Matrix3 newmat;
            newmat.set(1,1,get(2,2)*get(3,3) - get(2,3)*get(3,2));
            newmat.set(1,2,get(1,3)*get(3,2) - get(1,2)*get(3,3));
            newmat.set(1,3,get(1,2)*get(2,3) - get(1,3)*get(2,2));
            newmat.set(2,1,get(2,3)*get(3,1) - get(2,1)*get(3,3));
            newmat.set(2,2,get(1,1)*get(3,3) - get(1,3)*get(3,1));
            newmat.set(2,3,get(1,3)*get(2,1) - get(1,1)*get(2,3));
            newmat.set(3,1,get(2,1)*get(3,2) - get(2,2)*get(3,1));
            newmat.set(3,2,get(1,2)*get(3,1) - get(1,1)*get(3,2));
            newmat.set(3,3,get(1,1)*get(2,2) - get(1,2)*get(2,1));
            newmat /= det();
            return newmat;
        }

    };
}

#endif

