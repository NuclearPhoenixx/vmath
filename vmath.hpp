/*
    vmath.hpp

    Simple to use vector and matrix math library.

    2021, Phoenix1747, MIT License.
    
    Github: https://github.com/Phoenix1747/vmath

    TODO:
        - Orthonormalization -> Gram-Schmidt
        - Outer Product
        - Trace
*/

#ifndef VMATH_HPP
#define VMATH_HPP

#include <cmath> // Used only for pow, abs, sqrt, sin, cos, acos, M_PI
#include <vector>


namespace vmath{

    /*
    STANDARD N-DIMENSIONAL VECTOR FOR VECTOR MATH AND BASE FOR VECTOR2 & VECTOR3
    */
    template <class type=double> class VectorN{

        protected:
            std::vector<type> vec;

        public:
            VectorN(const std::size_t n=0){ // Constructor for vector dimension n
                vec.assign(n,0);
            }
            VectorN(const std::vector<type> &init){
                vec = init;
            }

            VectorN &operator=(const VectorN &newvec){
                if(this == &newvec) return *this;
                vec = newvec.vec;
                return *this; 
            }
            bool operator==(const VectorN &vector) const{
                return vec == vector.vec;
            }
            bool operator!=(const VectorN &vector) const{
                return !(vec == vector.vec);
            }
            bool operator<(const VectorN &vector) const{
                return length_squared() < vector.length_squared();
            }
            bool operator>(const VectorN &vector) const{
                return length_squared() > vector.length_squared();
            }
            bool operator<=(const VectorN &vector) const{
                return length_squared() <= vector.length_squared();
            }
            bool operator>=(const VectorN &vector) const{
                return length_squared() >= vector.length_squared();
            }
            type &operator[](const std::size_t pos){
                return vec[pos];
            }

            VectorN &operator+=(const VectorN &addvec){ // Vector Addition
                for(std::size_t i=0; i<vec.size(); i++) vec.at(i) += addvec.get(i);
                return *this;
            }
            VectorN &operator-=(const VectorN &addvec){ // Vector Subtraction
                for(std::size_t i=0; i<vec.size(); i++) vec.at(i) -= addvec.get(i);
                return *this;
            }
            VectorN &operator*=(const type num){ // Vector Scalar Multiplication
                for(std::size_t i=0; i<vec.size(); i++) vec.at(i) *= num;
                return *this;
            }
            VectorN &operator/=(const type num){ // Vector Scalar Division
                for(std::size_t i=0; i<vec.size(); i++) vec.at(i) /= num;
                return *this;
            }

            VectorN operator+(const VectorN &addvec) const{
                VectorN newvec = *this;
                newvec += addvec;
                return newvec;
            }
            VectorN operator-(const VectorN &addvec) const{
                VectorN newvec = *this;
                newvec -= addvec;
                return newvec;
            }
            VectorN operator*(const type num) const{
                VectorN newvec = *this;
                newvec *= num;
                return newvec;
            }
            VectorN operator/(const type num) const{
                VectorN newvec = *this;
                newvec /= num;
                return newvec;
            }

            std::size_t dim() const{ return vec.size(); } // Return vector dimension
            type get(const std::size_t i) const{ return type(vec.at(i)); }
            void set(const std::size_t i, type val){ vec.at(i) = val; }

            type length_squared() const{
                type l = 0;
                for(type num:vec) l += num*num;
                return l;
            }
            double length() const{ return std::sqrt(length_squared()); }
            
            bool is_normalized() const{
                return length_squared() == 1;
            }
            VectorN normalized() const{
                if(is_normalized()) return *this;
                VectorN newvec = *this;
                newvec /= newvec.length();
                return newvec;
            }
            type dot(const VectorN &dotvec) const{
                type dotprod = 0;
                for(std::size_t i=0; i<vec.size(); i++) dotprod += vec.at(i) * dotvec.get(i);
                return dotprod;
            }
            VectorN abs() const{
                std::vector<type> newvec;
                for(auto val:vec) newvec.push_back(std::abs(val));
                return VectorN(newvec);
            }
            type distance_squared_to(const VectorN &tovec) const{ return (*this-tovec).length_squared(); }
            double distance_to(const VectorN &tovec) const{ return (*this-tovec).length(); }
    };


    /*
    STANDARD 2D VECTOR FOR VECTOR MATH
    */
    template <class type=double> struct Vector2: public VectorN<type>{

        Vector2(const type initx=0, const type inity=0){ // Create vector with two coords
            this->vec = {initx, inity};
        }
        Vector2(const VectorN<type> &v){ // Create vector out of existing VectorN
            this->vec.assign(2,0);
            for(unsigned char i=0; i<2; i++) this->vec.at(i) = v.get(i);
        }

        type x() const{ return this->get(0); } // Return first coord
        type y() const{ return this->get(1); } // Return second coord

        double angle_to(const Vector2 &tovec) const{
            const double dotprod = this->dot(tovec);
            const double lengthprod = this->length() * tovec.length();
            return std::acos(dotprod/lengthprod);
        }
        double aspect() const{ return double(x()) / y(); }

        Vector2<double> rotated(const double phi) const{
            const double xi = x() * std::cos(phi) - y() * std::sin(phi);
            const double yi = y() * std::cos(phi) + x() * std::sin(phi);
            return Vector2<double>(xi,yi);
        }
        Vector2<double> tangent() const{ return rotated(M_PI/2); }
        
    };


    /*
    STANDARD 3D VECTOR FOR VECTOR MATH
    */
    template <class type=double> struct Vector3: public VectorN<type>{

        Vector3(const type initx=0, const type inity=0, const type initz=0){ // Create vector with the 3 coords
            this->vec = {initx, inity, initz};
        }
        Vector3(const VectorN<type> &v){ // Create vector out of existing vector
            this->vec.assign(3,0);
            for(unsigned char i=0; i<3; i++) this->vec.at(i) = v.get(i);
        }

        type x() const{ return this->get(0); }
        type y() const{ return this->get(1); }
        type z() const{ return this->get(2); }

        Vector3<double> cross(Vector3 &crossvec) const{
            const double newx = y() * crossvec.z() - z() * crossvec.y();
            const double newy = z() * crossvec.x() - x() * crossvec.z();
            const double newz = x() * crossvec.y() - y() * crossvec.x();
            return Vector3<double>(newx,newy,newz);
        }

    };


    /*
    STANDARD N-DIMENSIONAL MATRIX FOR VECTOR/MATRIX MATH AND BASE FOR MATRIX2 & MATRIX3
    */
    template <class type=double> class MatrixN{
        
        private:
            void remove_line(const std::size_t i=0){ mat.erase(mat.begin()+i); }
            void remove_column(const std::size_t i=0){ for(auto &matrix:mat) matrix.erase(matrix.begin()+i); }

            type subdet(const std::size_t i, const std::size_t j, MatrixN matrix) const{
                matrix.remove_line(i);
                matrix.remove_column(j);
                return matrix.det();
            }

        protected:
            std::vector<std::vector<type>> mat;

        public:
            MatrixN(std::size_t n=0, std::size_t m=0){ // Create (n x m) matrix with values 0
                mat.assign(n,std::vector<type>(m,0));
                while(n>0 && m>0){
                    mat.at(n-1).at(m-1) = 1;
                    n--;
                    m--;
                }
            }
            MatrixN(const std::vector<std::vector<type>> &init){ // Create matrix out of existing 2-dimensional vector
                mat = init;
            }
            
            MatrixN &operator=(const MatrixN &newmat){
                if(this == &newmat) return *this;
                mat.assign(newmat.dim().x(), std::vector<type>(newmat.dim().y(), 0));
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat.at(n).size(); m++){
                        mat.at(n).at(m) = newmat.get(n,m);
                    }
                }
                return *this; 
            }
            bool operator==(const MatrixN &matrix) const{
                return mat == matrix.mat;
            }
            bool operator!=(const MatrixN &matrix) const{
                return !(mat == matrix.mat);
            }
            std::vector<type> &operator[](const std::size_t pos){
                return mat[pos];
            }
            
            MatrixN &operator+=(const MatrixN &addmat){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat.at(n).size(); m++){
                        mat.at(n).at(m) += addmat.get(n, m);
                    }
                }
                return *this;
            }
            MatrixN &operator-=(const MatrixN &addmat){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat.at(n).size(); m++){
                        mat.at(n).at(m) -= addmat.get(n, m);
                    }
                }
                return *this;
            }
            MatrixN &operator*=(const type num){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat.at(n).size(); m++){
                        mat.at(n).at(m) *= num;
                    }
                }
                return *this;
            }
            MatrixN &operator*=(const MatrixN &multmat){
                MatrixN tempmat(dim().x(), dim().y());
                type val = 0;

                for(std::size_t i=0; i<mat.size(); i++){
                    for(std::size_t k=0; k<mat.at(i).size(); k++){
                        val = 0;
                        for(std::size_t j=0; j<mat.at(i).size(); j++){
                            val += mat.at(i).at(j) * multmat.get(j,k);
                        }
                        tempmat.set(i,k,val);
                    }
                }
                *this = tempmat; // meh a bit ugly
                return *this;
            }
            MatrixN &operator/=(const type num){
                for(std::size_t n=0; n<mat.size(); n++){
                    for(std::size_t m=0; m<mat.at(n).size(); m++){
                        mat.at(n).at(m) /= num;
                    }
                }
                return *this;
            }
            
            MatrixN operator+(const MatrixN &addmat) const{
                MatrixN newmat = *this;
                newmat += addmat;
                return newmat;
            }
            MatrixN operator-(const MatrixN &addmat) const{
                MatrixN newmat = *this;
                newmat -= addmat;
                return newmat;
            }
            MatrixN operator*(const type num) const{ // Matrix Scalar Operation
                MatrixN newmat = *this;
                newmat *= num;
                return newmat;
            }
            MatrixN operator*(const MatrixN &multmat) const{ // Matrix Matrix Operation
                MatrixN newmat = *this;
                newmat *= multmat;
                return newmat;
            }
            VectorN<type> operator*(const VectorN<type> &multvec) const{ // Matrix Vector Operation
                VectorN<type> newvec(multvec.dim());
                type val = 0;

                for(std::size_t i=0; i<multvec.dim(); i++){
                    val = 0;
                    for(std::size_t j=0; j<dim().y(); j++){
                        val += get(i,j) * multvec.get(j);
                    }
                    newvec.set(i,val);
                }
                return newvec;
            }
            MatrixN operator/(const type num) const{
                MatrixN newmat = *this;
                newmat /= num;
                return newmat;
            }

            Vector2<std::size_t> dim() const{
                if(mat.size() == 0) return Vector2<std::size_t>();
                return Vector2<std::size_t>(mat.size(),mat.at(0).size());
            }
            type get(const std::size_t n, const std::size_t m) const{ return mat.at(n).at(m); }
            void set(const std::size_t n, const std::size_t m, type val) { mat.at(n).at(m) = val; }

            type det(const std::size_t i=1, const std::size_t j=1) const{ // Todo: Change starting index back to 0.
                if(dim().x() == 2) return get(0,0) * get(1,1) - get(1,0) * get(0,1);

                type val = 0;
                for(std::size_t jj=j; jj<=dim().y(); jj++){
                    val += std::pow(-1,i+jj) * get(i-1,jj-1) * subdet(i-1,jj-1,*this);
                }
                return val;
            }
            MatrixN transpose() const{
                MatrixN<type> newmat = *this;
                for(std::size_t n=0; n<dim().x(); n++){
                    for(std::size_t m=0; m<dim().y(); m++) newmat.set(n,m,get(m,n));
                }
                return newmat;
            }
            bool is_normalized() const{ // If length of every column == 1
                for(std::size_t m=0; m<dim().y(); m++){
                    std::vector<type> tempvec;
                    for(std::size_t n=0; n<dim().x(); n++) tempvec.push_back(get(n,m));
                    if(!VectorN<type>(tempvec).is_normalized()) return false;
                }
                return true;
            }
            bool is_orthogonalized() const{ // If every column is perpendicular to every other one
                std::vector<VectorN<type>> vecholder;
                
                for(std::size_t m=0; m<dim().y(); m++){
                    std::vector<type> tempvec;
                    for(std::size_t n=0; n<dim().x(); n++) tempvec.push_back(get(n,m));
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
            bool is_orthonormalized() const{
                return is_orthogonalized() && is_normalized();
            }
            
    };


    /*
    STANDARD 2x2 MATRIX FOR VECTOR MATH
    */
    template <class type=double> struct Matrix2: public MatrixN<type>{

        Matrix2(type a11=1, type a12=0, type a21=0, type a22=1){ // Construct with matrix values
            this->mat = {{a11,a12},{a21,a22}};
        }
        Matrix2(const std::vector<std::vector<type>> &init){ // Create matrix out of existing 2x2 vector
            this->mat = init;
        }
        Matrix2(const MatrixN<type> &init){ // Construct with existing 2x2 MatrixN
            this->mat.assign(2, std::vector<type>(2,0));
            for(unsigned char n=0; n<2; n++){
                for(unsigned char m=0; m<2; m++){
                    this->mat.at(n).at(m) = init.get(n,m);
                }
            }
        }
        
        Matrix2 inverse() const{
            Matrix2 newmat;
            newmat.set(0,0, this->get(1,1));
            newmat.set(0,1, -this->get(0,1));
            newmat.set(1,0, -this->get(1,0));
            newmat.set(1,1, this->get(0,0));
            newmat /= this->det();
            return newmat;
        }

    };


    /*
    STANDARD 3x3 MATRIX FOR VECTOR MATH
    */
    template <class type=double> struct Matrix3: public MatrixN<type>{

        Matrix3(type a11=1, type a12=0, type a13=0, type a21=0, type a22=1, type a23=0, type a31=0, type a32=0, type a33=1){ // Construct with matrix values
            this->mat = {{a11,a12,a13},{a21,a22,a23},{a31,a32,a33}};
        }
        Matrix3(const std::vector<std::vector<type>> &init){ // Construct with 2d 3x3 vector
            this->mat = init;
        }
        Matrix3(const MatrixN<type> &init){ // Construct with existing 3x3 MatrixN
            this->mat.assign(3, std::vector<type>(3,0));
            for(unsigned char n=0; n<3; n++){
                for(unsigned char m=0; m<3; m++){
                    this->mat.at(n).at(m) = init.get(n,m);
                }
            }
        }
        
        Matrix3 inverse() const{
            Matrix3 newmat;
            newmat.set(0,0, this->get(1,1) * this->get(2,2) - this->get(1,2) * this->get(2,1));
            newmat.set(0,1, this->get(0,2) * this->get(2,1) - this->get(0,1) * this->get(2,2));
            newmat.set(0,2, this->get(0,1) * this->get(1,2) - this->get(0,2) * this->get(1,1));
            newmat.set(1,0, this->get(1,2) * this->get(2,0) - this->get(1,0) * this->get(2,2));
            newmat.set(1,1, this->get(0,0) * this->get(2,2) - this->get(0,2) * this->get(2,0));
            newmat.set(1,2, this->get(0,2) * this->get(1,0) - this->get(0,0) * this->get(1,2));
            newmat.set(2,0, this->get(1,0) * this->get(2,1) - this->get(1,1) * this->get(2,0));
            newmat.set(2,1, this->get(0,1) * this->get(2,0) - this->get(0,0) * this->get(2,1));
            newmat.set(2,2, this->get(0,0) * this->get(1,1) - this->get(0,1) * this->get(1,0));
            newmat /= this->det();
            return newmat;
        }

    };

}

#endif

