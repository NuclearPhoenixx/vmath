# vmath.hpp

C++ Library introducing Vector2, Vector3, Matrix2 and Matrix3 objects as well as n-dimensional VectorN and MatrixN to be used in vector/matrix math.

For documentation please head to the [docs sheet](DOCS.md). All classes are inside the namespace `vmath`.

Needs following parts of the standard library:
- std::cmath for a couple of math functions and the Pi constant.
- std::vector for the underlying data structure.

## Features
- 2- and 3-dimensional as well as arbitrary n-dimensional objects.
- Reasonable matrix and vector comparison and assignment operators.
- Reasonable arithmetic operations using operators.
- Basic linear algebra operations like dot product, cross product, transpose (only matrices), inverse (only Matrix2 and Matrix3), determinant...

