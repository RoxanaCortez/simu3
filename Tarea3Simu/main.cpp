#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"
#include "display_tools.h"

int main(){

    vector<Matrix> localKs;
    vector<Vector> localbs;
    // ahora nuestra incognita es el vector U
    Matrix K;
    Vector b;
    Vector U;

    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- 1 DIMENSI"<<char(224)<<"N\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m);

    crearSistemasLocales(m,localKs,localbs);

    zeroes(K,m.getSize(NODES));
    zeroes(b,m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);

    applyNeumann(m,b);

    applyDirichlet(m,K,b);
    zeroes(U,b.size());
    calculate(K,b,U);

    cout << "La respuesta es: \n";
    showVector(U);

    return 0;
}
