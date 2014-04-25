#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592


FLOAT * get_memory(int n_points){
    FLOAT * x;
    if(!(x = malloc(sizeof(FLOAT) * n_points))){
        printf("problem with memory allocation");
        exit(1);
    }
    return x;
}

void copy(FLOAT *u, FLOAT *u_past, int n_points){
    
    int i;
    for(i=0;i<n_points;i++){
        u_past[i] = u[i];
    }
}



int main(int argc, char **argv){

    
    int n_points = 1000; //n_points in x
    int n_t = 350;
    
    FLOAT *x;
    FLOAT *u;
    FLOAT *u_past;
    FLOAT *u_initial;
    FLOAT *u_future;
    FLOAT *u_present;
    FLOAT r;
    FLOAT dx;
    FLOAT dt;
    FLOAT c = 1.0;
    x = get_memory(n_points);
    u = get_memory(n_points);
    u_initial = get_memory(n_points);
    u_past = get_memory(n_points);
    u_future = get_memory(n_points);
    u_present = get_memory(n_points);
    int i,j;
    
    x[0] = 0.0;
    
    for (i=0; i<n_points; i++) {
        
        if (i!=0) {
        x[i] = x[i-1] + (1.0)/n_points;
        }
        
        u_initial[i] = exp(-((x[i]-0.3)*(x[i]-0.3))/0.01);
    }
    
    
    /*find the first iteration*/
    
    dx = x[1]-x[0];
    dt = 0.0005;
    r = c * dt/dx;
    
    /*r should be less than 1.0 in this scheme in order to be numerically stable*/
    printf("%f\n",r);
    
    /*Fixed boundary condition*/
    u_initial[0] = 0.0;
    u_initial[n_points-1] = 0.0;
    
    u_future[0] = 0.0;
    u_future[n_points - 1] = 0.0;
    
    
    for (i=1; i<n_points - 1; i++) {
        u_future[i] = u_initial[i] + (pow(r,2)/2.0)*(u_initial[i+1]-2.0*u_initial[i]+ u_initial[i-1]);
    }
    //Variable to hold the previous value
    copy(u_initial, u_past, n_points);
    //Variable to hold the present value
    copy(u_future, u_present, n_points);
    
    
    /*Next iterations:*/
    
    for (j=0; j<n_t; j++) {
        printf("Voy en el paso %d \n",j);
        for (i=1; i<n_points; i++) {
            u_future[i] = (2.0*(1.0-pow(r,2)))*u_present[i] - u_past[i] + pow(r,2)*(u_present[i+1] + u_present[i-1]);
            printf("%f\n",u_present[i]);
        }
        //Variable to hold the previous value
        copy(u_present, u_past, n_points);
        //Variable to hold the present value
        copy(u_future, u_present, n_points);
    }
    
    
    FILE *in = fopen("datos.dat","w");
    
    
    for(i=0;i<n_points;i++){
        fprintf(in,"%f %f %f %f \n",x[i],u_initial[i], u_future[i],u_present[i]);
    }
    
    fclose(in);
    
    return 0;
    
    
}
