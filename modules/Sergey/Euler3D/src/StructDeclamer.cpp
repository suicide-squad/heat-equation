//
// Created by lenferd on 14.03.17.
//

#include "StructDeclamer.h"

int setTaskExpr(TaskExpressions &taskExpr, Task2 &task2){
    taskExpr.x1 = (task2.sigma * task2.dt) / (task2.timeStepX * task2.timeStepX);
//    printf("sigma %lf\n", sigma);
//    printf("dt %lf\n", dt);
//    printf("timeStepX %lf\n", timeStepX);
//    printf("taskx1 %lf\n", task.x1);
    taskExpr.y1 = (task2.sigma * task2.dt) / (task2.timeStepY * task2.timeStepY);
//    printf("tasky1 %lf\n", task.y1);
    taskExpr.z1 = (task2.sigma * task2.dt) / (task2.timeStepZ * task2.timeStepZ);
//    printf("taskz1 %lf\n", task.z1);

    taskExpr.x2Comp = (1 - 2 * taskExpr.x1 - 2 * taskExpr.y1 - 2 * taskExpr.z1);
//    printf("taskx2 %lf\n", task.x2Comp);

//1 - 2 * ((sigma * dt) * 1 / (timeStepX * timeStepX) * 1 / (timeStepY * timeStepY) * 1/ (timeStepZ * timeStepZ);

}