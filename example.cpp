/*
    Example testfile to experiment with vmath.hpp
*/

#include "vmath.hpp"
#include <iostream>

using namespace vmath;

int main(){
    // Basic VectorN, Vector2, Vector3 Types    
    VectorN<> va({1,2});
    VectorN<> vb(va);
    Vector2<> vc(2,-2);
    Vector3<> ZERO;

	// Simple Vector Operations
    std::cout << "ZERO.dim(): " << ZERO.dim() << std::endl;
    std::cout << "va[1]: " << va[1] << std::endl;
    std::cout << "(va > vb): " << (va > vb) << std::endl;
	
    // Vector Matrix Operations
    Vector2<> vd(vb);
    Matrix2<> mata(1,2,3,4);
    std::cout << "(mata * vd).dot(vc): " << (mata * vd).dot(vc) << std::endl;
    
	// Vector Comparison
	std::cout << "va == vb: " << (va==vb) << std::endl;
    
	// More Vector Operations
    std::cout << "vd + vb: " << Vector2<>(vd+vb).x() << "," << Vector2<>(vd+vb).y() << std::endl;    
    std::cout << "vd * 4: " << Vector2<>(vd*4).x() << "," << Vector2<>(vd*4).y() << std::endl;
    std::cout << "vd/2: " << Vector2<>(vd/2).x() << "," << Vector2<>(vd/2).y() << std::endl;
    std::cout << "(vd*4).normalized(): " << Vector2<>((vd*4).normalized()).x() << "," << Vector2<>((vd*4).normalized()).y() << std::endl;
    std::cout << "(vd/2).length(): " << ((vd/2).length()) << std::endl;
    std::cout << "vd: " << vd.x() << "," << vd.y() << std::endl;
    std::cout << "Vector2 vb: " << Vector2<>(vb).x() << "," << Vector2<>(vb).y() << std::endl;
    std::cout << "vd.length(): " << vd.length() << std::endl;
    std::cout << "vd.is_normalized(): " << vd.is_normalized() << std::endl;
    
	// Vector3 Cross Product
    Vector3<unsigned int> aa(1,2,3);
    Vector3<unsigned int> bb(4,6,2);
    Vector3<> dd = aa.cross(bb);

    std::cout << "aa.cross(bb): " << dd.x() << "," << dd.y() << "," << dd.z() << std::endl;
        
    return 0;
}

