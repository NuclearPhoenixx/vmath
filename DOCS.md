# class VectorN

```cpp
std::vector<double> vec;
```
Underlying datatype everything is built upon. Holds all vector values.

```cpp
VectorN(double n=0.)
```
Constructor for n-dimensional Vector with values 0.

```cpp
VectorN(const std::vector<double> &init)
```
Constructor taking all the values of an existing vector.

```cpp
unsigned int dim()
```
Returns the vector dimension.

```cpp 
double get(unsigned int i)
```
Returns value at position i. Starting at 1.

```cpp
void set(unsigned int i, double val)
```
Sets vector value at positon i to val.

```cpp
double length_squared()
```
Returns length of the vector squared.

```cpp
double length()
```
Returns length of the vector.

```cpp
bool is_normalized()
```
Returns true when the vector length is equal to 1.

```cpp
VectorN normalized()
```
Returns another vector with length equal to 1.

```cpp
double dot(const VectorN &dotvec)
```
Returns the dot product of the vector and dotvec.

```cpp
VectorN abs()
```
Returns another vector with all coordinates being positive.

```cpp
double distance_squared_to(const VectorN &tovec)
```
Returns the distance between the vector and tovec squared.

```cpp
double distance_to(const VectorN &tovec)
```
Returns the distance between the vector and tovec.

# struct Vector2

**Vector2 inherits all functionality from VectorN. Therefore, you can do stuff like `Vector2.length()` just like above.**

**Vector2 features following additional properties and functions:**

```cpp
VectorN vec;
```
Vector2 is based on this generic VectorN.

```cpp
Vector2(const double initx=0., const double inity=0.)
```
Constructor with given values.

```cpp
Vector2(const VectorN &v)
```
Construct a Vector2 out of a _2-dimensional_ VectorN.

```cpp
double x()
```
Returns the x (first) coordinate.

```cpp
double y()
```
Returns the y (second) coordinate.

```cpp
double angle_to(const Vector2 &tovec)
```
Returns the angle between the vector and tovec in radians.

```cpp
double aspect()
```
Returns the ratio of x to y coordinate. Y must be non-zero!

```cpp
Vector2 rotated(double phi)
```
Returns another vector rotated counter-clockwise by phi radians.

```cpp
Vector2 tangent()
```
Returns another vector that is perpendicular to the vector (aka rotated by 90°).

# struct Vector3

**Vector3 inherits all functionality from VectorN. Therefore, you can do stuff like `Vector3.length()` just like above.**

**Vector3 features following additional properties and functions:**

```cpp
VectorN vec;
```
Vector3 is based on this generic VectorN.

```cpp
Vector3(const double initx=0., const double inity=0., const double initz=0.)
```
Constructs a new Vector3 with values initx, inity and initz.

```cpp
Vector3(const VectorN &v)
```
Constructs a Vector3 out of a _3-dimensional_ VectorN.

```cpp
double x()
```
Returns the x (first) coordinate.

```cpp
double y()
```
Returns the y (second) coordinate.

```cpp
double z()
```
Returns the z (third) coordinate.

```cpp
Vector3 cross(Vector3 &crossvec)
```
Returns a new vector that is perpendicular to the vector and crossvec.

# class MatrixN

```cpp
MatrixN(double n=0.,double m=0.)
```
Construct a (n x m) matrix with values all being 0.

```cpp
MatrixN(const std::vector<std::vector<double>> &init)
```
Constructs a matrix out of an existing two-dimensional vector.

```cpp
Vector2 dim()
```
Returns a Vector2 with the matrix dimensions. x is n dimension, y is m dimension of a (n x m) matrix.

```cpp
double get(unsigned int n, unsigned int m)
```
Gets the value at the nth row and mth column. 

```cpp
void set(unsigned int n, unsigned int m, double val)
```
Sets the value at the nth row and mth column to val.

```cpp
double det(unsigned int i=1, unsigned int j=1)
```
Returns the determinant. Uses Laplace expansion starting at n=i=1, m=j=1 default.

```cpp
MatrixN transpose()
```
Returns a new matrix that is the transpose of the matrix.

```cpp
bool is_normalized()
```
Returns true if each column of the matrix has a length of 1 when seen as a VectorN.

```cpp
bool is_orthogonalized()
```
Returns true if each column of the matrix is perpendicular to every other column. 

```cpp
bool is_orthonormalized()
```
Returns true if the matrix is normalized and orthogonalized.

# struct Matrix2

**Matrix2 inherits all functionality from MatrixN. Therefore, you can do stuff like `Matrix2.transpose()` just like above.**

**Matrix2 features following additional properties and functions:**

```cpp
MatrixN mat;
```
Matrix2 is built upon a generic MatrixN object. This holds all the values.

```cpp
Matrix2(double a11=0., double a12=0., double a21=0., double a22=0.)
```
Construct a Matrix2 object with all the according values.

```cpp
Matrix2(const std::vector<std::vector<double>> &init)
```
Construct a Matrix2 object from a 2x2 vector.

```cpp
Matrix2(const MatrixN &init)
```
Construct a Matrix2 object from a 2x2 MatrixN object.

```cpp
Matrix2 inverse()
```
Returns a new matrix which is the inverse of the matrix.

# struct Matrix3

**Matrix3 inherits all functionality from MatrixN. Therefore, you can do stuff like `Matrix3.transpose()` just like above.**

**Matrix3 features following additional properties and functions:**

```cpp
MatrixN mat;
```
Matrix3 is built upon a generic MatrixN object. This holds all the values.

```cpp
Matrix3(double a11=0.,double a12=0.,double a13=0.,double a21=0.,double a22=0.,double a23=0.,double a31=0.,double a32=0.,double a33=0.)
```
Construct a Matrix3 object with all the according values.

```cpp
Matrix3(const std::vector<std::vector<double>> &init)
```
Construct a Matrix3 object from a 3x3 vector.

```cpp
Matrix3(const MatrixN &init)
```
Construct a Matrix3 object from a 3x3 MatrixN object.

```cpp    
Matrix3 inverse()
```
Returns a new matrix which is the inverse of the matrix.

