# VectorN
### `template <class type=double> class VectorN`

```cpp
VectorN(const std::size_t n=0)
```
Constructor for n-dimensional Vector with values 0.

```cpp
VectorN(const std::vector<type> &init)
```
Constructor taking all the values of an existing vector.

```cpp
VectorN(const VectorN<type> &init)
```
Copy Constructor

```cpp
std::size_t dim()
```
Returns the vector dimension.

```cpp 
type get(const std::size_t i)
```
Returns value at position i. Starting at i=0. With bounds checking.

```cpp
void set(const std::size_t i, type val)
```
Sets vector value at positon i to val. With bounds checking.

```cpp
type length_squared()
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
type dot(const VectorN &dotvec)
```
Returns the dot product of the vector and dotvec.

```cpp
VectorN abs()
```
Returns another vector with all coordinates being positive.

```cpp
type distance_squared_to(const VectorN &tovec)
```
Returns the distance between the vector and tovec squared.

```cpp
double distance_to(const VectorN &tovec)
```
Returns the distance between the vector and tovec.

### Comparison Operators

```cpp
bool operator==(const VectorN &vector)
bool operator!=(const VectorN &vector)
```
Checks if two vectors are _exactly_ the same vector.

```cpp
bool operator<(const VectorN &vector)
bool operator>(const VectorN &vector)
bool operator<=(const VectorN &vector)
bool operator>=(const VectorN &vector)
```
Compares the length of two vectors (uses length_squared() for performance).

```cpp
type &operator[](const std::size_t pos)
```
Returns a reference to the element at pos. No bounds checking.

# Vector2
### `template <class type=double> struct Vector2: public VectorN<type>`

**Vector2 inherits all functionality from VectorN. Therefore, you can do stuff like `Vector2.length()` just like above.**

**Vector2 features following additional properties and functions:**

```cpp
Vector2(const type initx=0, const type inity=0)
```
Constructor with given values.

```cpp
Vector2(const VectorN<type> &v)
```
Construct a Vector2 out of a _2-dimensional_ VectorN.

```cpp
type x()
```
Returns the x (first) coordinate.

```cpp
type y()
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
Vector2<double> rotated(const double phi)
```
Returns another vector (always double!) rotated counter-clockwise by phi radians.

```cpp
Vector2<double> tangent()
```
Returns another vector that is perpendicular to the vector. This is an alias for rotated(M_PI/2).

# Vector3
### `template <class type=double> struct Vector3: public VectorN<type>`

**Vector3 inherits all functionality from VectorN. Therefore, you can do stuff like `Vector3.length()` just like above.**

**Vector3 features following additional properties and functions:**

```cpp
Vector3(const type initx=0, const type inity=0, const type initz=0)
```
Constructs a new Vector3 with values initx, inity and initz.

```cpp
Vector3(const VectorN<type> &v)
```
Constructs a Vector3 out of a _3-dimensional_ VectorN.

```cpp
type x()
```
Returns the x (first) coordinate.

```cpp
type y()
```
Returns the y (second) coordinate.

```cpp
type z()
```
Returns the z (third) coordinate.

```cpp
Vector3<double> cross(Vector3 &crossvec)
```
Returns a new vector (always double!) that is perpendicular to the vector and crossvec.

# MatrixN
### `template <class type=double> class MatrixN`

```cpp
MatrixN(const std::size_t n=0, const std::size_t m=0)
```
Construct an (n x m) matrix. Defaults to the identity matrix if dimensions are non-zero.

```cpp
MatrixN(const std::vector<std::vector<type>> &init)
```
Constructs a matrix out of an existing two-dimensional vector.

```cpp
Vector2<std::size_t> dim()
```
Returns a Vector2 (always std::size_t!) with the matrix dimensions. x is n dimension, y is m dimension of a (n x m) matrix.

```cpp
type get(const std::size_t n, const std::size_t m)
```
Gets the value at the nth row and mth column. Starts at 0. With bounds checking.

```cpp
void set(const std::size_t n, const std::size_t m, type val)
```
Sets the value at the nth row and mth column to val. With bounds checking.

```cpp
type det(const std::size_t i=1, const std::size_t j=1)
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

### Comparison Operators

```cpp
bool operator==(const MatrixN &matrix)
bool operator!=(const MatrixN &matrix)
```
Checks if two matrices are _exactly_ the same matrix.

```cpp
std::vector<type> &operator[](const std::size_t pos)
```
Returns an std::vector reference to the row at pos. No bounds checking.

# Matrix2
### `template <class type=double> struct Matrix2: public MatrixN<type>`

**Matrix2 inherits all functionality from MatrixN. Therefore, you can do stuff like `Matrix2.transpose()` just like above.**

**Matrix2 features following additional properties and functions:**

```cpp
Matrix2(type a11=1, type a12=0, type a21=0, type a22=1)
```
Construct a Matrix2 object with all the according values. Defaults to the identity matrix.

```cpp
Matrix2(const std::vector<std::vector<type>> &init)
```
Construct a Matrix2 object from a 2x2 vector.

```cpp
Matrix2(const MatrixN<type> &init)
```
Construct a Matrix2 object from a 2x2 MatrixN object.

```cpp
Matrix2 inverse()
```
Returns a new matrix which is the inverse of the matrix.

# Matrix3
### `template <class type=double> struct Matrix3: public MatrixN<type>`

**Matrix3 inherits all functionality from MatrixN. Therefore, you can do stuff like `Matrix3.transpose()` just like above.**

**Matrix3 features following additional properties and functions:**

```cpp
Matrix3(type a11=1, type a12=0, type a13=0, type a21=0, type a22=1, type a23=0, type a31=0, type a32=0, type a33=1)
```
Construct a Matrix3 object with all the according values. Defaults to the identity matrix.

```cpp
Matrix3(const std::vector<std::vector<type>> &init)
```
Construct a Matrix3 object from a 3x3 vector.

```cpp
Matrix3(const MatrixN<type> &init)
```
Construct a Matrix3 object from a 3x3 MatrixN object.

```cpp    
Matrix3 inverse()
```
Returns a new matrix which is the inverse of the matrix.

