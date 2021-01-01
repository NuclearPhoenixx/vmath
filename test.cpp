/*
    Testfile to experiment with vmath.hpp
*/

#include "vmath.hpp"
#include <iostream>
#include <cmath>
#include <vector>

using namespace vmath;
//using namespace std;

int main(){
    Vector2 z(-2,3);
    Vector2 y(1,0);
    std::cout << y.tangent().x() << "," << y.tangent().y() << std::endl;
    /*
    Vector2 i(1,0);
    Vector2 a(1,1);
    Vector2 b(2,2);
    Vector2 d = a;
    
    std::cout << (a+b).x() << "," << (a+b).y() << std::endl;
    std::cout << (a*4).x() << "," << (a*4).y() << std::endl;
    std::cout << (a/2).x() << "," << (a/2).y() << std::endl;
    std::cout << ((a*4).normalized()).x() << "," << ((a*4).normalized()).y() << std::endl;
    std::cout << ((a/2).length()) << std::endl;
    std::cout << a.x() << "," << a.y() << std::endl;
    std::cout << b.x() << "," << b.y() << std::endl;
    std::cout << ((a).length()) << std::endl;
    
    d += Vector2(2,2);
    std::cout << d.x() << "," << d.y() << std::endl;
    
    d *= 3;
    std::cout << d.x() << "," << d.y() << std::endl;
    d /= 3;
    std::cout << d.x() << "," << d.y() << std::endl;
    d -= a;
    std::cout << d.x() << "," << d.y() << std::endl;
    std::cout << a.x() << "," << a.y() << std::endl;

    double dot = a.dot(b);
    std::cout << dot << std::endl;
    std::cout << i.is_normalized() << std::endl;
    
    Vector3 aa(1,2,3);
    Vector3 bb(4,6,2);
    Vector3 dd = aa.cross(bb);

    std::cout << dd.x() << "," << dd.y() << "," << dd.z() << std::endl;
        
    VectorN aN;
    VectorN bN(4);
    std::vector<double> vectorr = {1,2,3,4,5};
    VectorN cN(vectorr);

    std::cout << aN.dim() << std::endl;
    std::cout << bN.dim() << std::endl;
    std::cout << cN.dim() << std::endl;
    std::cout << cN.get(5) << std::endl;
    */    

    std::vector<std::vector<double>> lola = {{1,2},{3,4}};
    std::vector<std::vector<double>> lolb = {{5,6},{7,8}}; 
    MatrixN c(lola);
    MatrixN b(lolb);
    //a *= b;
    MatrixN a = b * c;
    /*
    Vector2 dim = c.dim();
    for(size_t n=1; n<=dim.x(); n++){
        for(size_t m=1; m<=dim.y(); m++){
            std::cout << c.get(n,m) << ",";
        }
        std::cout << std::endl;
    }
    */
    VectorN anotherTest({1,2,3});
    MatrixN test({{1,0,0},{0,1,0},{0,0,1}});
    //std::cout << test.det() << std::endl;
    
    VectorN thingy = test * (test+test) * anotherTest;

    for(size_t i=1; i<=thingy.dim(); i++) std::cout << thingy.get(i) << std::endl;
    /*
    Vector2 dim = test.dim();
    for(size_t n=1; n<=dim.x(); n++){
        for(size_t m=1; m<=dim.y(); m++){
            std::cout << test.get(n,m) << ",";
        }
        std::cout << std::endl;
    }*/
    /*dim = c.dim();
    for(size_t n=1; n<=dim.x(); n++){
        for(size_t m=1; m<=dim.y(); m++){
            std::cout << c.get(n,m) << ",";
        }
        std::cout << std::endl;
    }
    dim = b.dim();
    for(size_t n=1; n<=dim.x(); n++){
        for(size_t m=1; m<=dim.y(); m++){
            std::cout << b.get(n,m) << ",";
        }
        std::cout << std::endl;
    }
    */

    return 0;
}
