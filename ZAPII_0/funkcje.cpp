#include "funkcje.h"

void run()
{
    Matrix matA,matB,matX,l,u,c;

    matA.open("A.txt");
    matA.fill("A.txt");
    matB.open("B.txt");
    matB.fill("B.txt");
    matA.lu(l,u);
    l.invert();
    l.save("l.txt");
    u.invert();
    u.save("u.txt");
    c.multiple(l,matB);
    c.save("c.txt");
    matX.multiple(u,c);
    matX.save("SaveX.txt");


//    mat1.open("test.txt");
//    mat1.fill("test.txt");
//    mat2.open("savetest2.txt");
//    mat2.fill("savetest2.txt");
//    mat.multiple(mat1,mat2);
//    mat.save("saveTest3.txt");
//    mat1.transp();
//    mat1.save("saveTransp.txt");
//    mat1.lu(l,u);
//    l.save("saveL.txt");
//    u.save("saveU.txt");

}

