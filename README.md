# vmath.hpp

C++ Library introducing Vector2, Vector3, Matrix2 and Matrix3 objects as well as n-dimensional VectorN and MatrixN to be used in vector/matrix math.

To include this lib in your project you only need to grab the [vmath.hpp](vmath.hpp) file and you're ready to go!

* For documentation please head to the [docs sheet](DOCS.md).
* All classes are inside the namespace `vmath`.
* Only use Vector2 with other Vector2 or Matrix2 objects. The same is true for Vector3/Matrix3 and VectorN/MatrixN. In general **never mix object classes of different dimensions** - that's very important, obviously! **However**, you could e.g. mix a 2-D VectorN with a Vector2 or a Matrix2. To be clear, this is not the recommended use case and some things might break, but most of the basic functionality *should* work.

Needs following parts of the C++ standard library:
- `<cmath>` for a couple of math functions and the Pi constant.
- `<vector>` for the underlying data structure.

## Features
- 2- and 3-dimensional as well as arbitrary n-dimensional objects.
- Vector and matrix classes support different data types.
- Reasonable matrix and vector comparison and assignment operators.
- Reasonable arithmetic operations using operators.
- Basic linear algebra operations like dot product, cross product, transpose (only matrices), inverse (only Matrix2 and Matrix3), determinant...

