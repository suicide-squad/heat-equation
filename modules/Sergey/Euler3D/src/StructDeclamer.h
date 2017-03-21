//
// Created by lenferd on 14.03.17.
//

#ifndef EULER3D_STRUCTDECLAMER_H
#define EULER3D_STRUCTDECLAMER_H

struct TaskExpressions {
    double x1;
    double x2Comp;

    double y1;


    double z1;
//    double z2;
} ;

struct SparseMatrix {
    int _size;
    int _rows;
    double *values;
    int *columns;   // какой столбец
    int *pointerB;  // указатель на начало строки
};

#endif //EULER3D_STRUCTDECLAMER_H
