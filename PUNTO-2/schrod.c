#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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

    FILE *in;
    in = fopen("datos.dat","w");
    
    int n_points = 500; //n_points in x
    int n_t = atoi(argv[1]);
    
    
    /*Initializes all the variables that we will need*/
    
    
    FLOAT r;
    FLOAT dx;
    FLOAT dt;
    FLOAT L = 20.0;
    FLOAT k_0 = 16*PI;
    FLOAT sigma_0 = 0.05;
    FLOAT *x = get_memory(n_points);
    FLOAT *Re_u = get_memory(n_points);
    FLOAT *Re_u_initial = get_memory(n_points);
    FLOAT *Re_u_past = get_memory(n_points);
    FLOAT *Re_u_future = get_memory(n_points);
    FLOAT *Re_u_present = get_memory(n_points);
    
    FLOAT *Im_u = get_memory(n_points);
    FLOAT *Im_u_initial = get_memory(n_points);
    FLOAT *Im_u_past = get_memory(n_points);
    FLOAT *Im_u_future = get_memory(n_points);
    FLOAT *Im_u_present = get_memory(n_points);
    
    FLOAT *rho = get_memory(n_points);
    
    FLOAT *V = get_memory(n_points);
    int i,j;
    
    /*Initializes the conditions with a Gaussian*/
    
    x[0] = -5.0;
    
    for (i=0; i<n_points; i++) {
        
        if (i!=0) {
        x[i] = x[i-1] + L/n_points;
        }
    
        
        Re_u_initial[i] = exp(-(1.0/2.0)*pow((1/sigma_0)*(x[i]-5.0),2))*cos(k_0*x[i]);
        Im_u_initial[i] = exp(-(1.0/2.0)*pow((1/sigma_0)*(x[i]-5.0),2))*sin(k_0*x[i]);
        
        V[i] = 0.5*pow(x[i],2);
        
        
        fprintf(in,"%f %f %f \n",x[i], V[i], pow(Re_u_initial[i],2)+pow(Im_u_initial[i],2));
        
    }
    fclose(in);

    /*find the numerical stability and the first iteration*/
    
    dx = x[1]-x[0];
    dt = 0.001;
    FLOAT alpha = dt/(2*pow(dx,2));
    
    printf("%f %f dt dx \n", dt, dx);
    
    /*alpha should be less than 0.5 in this scheme in order to be numerically stable*/
    printf("Este valor es alpha: alpha vale %f y debe ser menor a 0.5 \n",alpha);

    FILE *in2 = fopen("datosevol.dat","w");
    
    
    /*Fixed boundary condition*/
    Re_u_initial[0] = 0.0;
    Re_u_initial[n_points-1] = 0.0;
    Im_u_initial[0] = 0.0;
    Im_u_initial[n_points-1] = 0.0;
    
    Re_u_future[0] = 0.0;
    Re_u_future[n_points - 1] = 0.0;
    Im_u_future[0] = 0.0;
    Im_u_future[n_points - 1] = 0.0;
    
    
    for (i=1; i<n_points - 1; i++) {
        Re_u_present[i] = Re_u_initial[i]-2*(alpha*(Im_u_initial[i+1]+Im_u_initial[i-1])-2*(alpha+V[i]*dt)*Im_u_initial[i]);
        Im_u_present[i] = Im_u_initial[i]+2*(alpha*(Re_u_initial[i+1]+Re_u_initial[i-1])-2*(alpha+V[i]*dt)*Re_u_initial[i]);
        
        
        rho[i] = pow(Re_u_present[i],2) + Im_u_present[i]*Im_u_initial[i];
        fprintf(in2,"%f ", rho[i]);
        
    }
    
    fprintf(in2, "\n");
    
    //Variables to hold the previous value
    copy(Re_u_initial, Re_u_past, n_points);
    copy(Im_u_initial, Im_u_past, n_points);
    //Variables to hold the present value
    copy(Re_u_future, Re_u_present, n_points);
    copy(Im_u_future, Im_u_present, n_points);
    
    
    
    
    /*Next iterations:*/
    
    for (j=1; j<n_t; j++) {
        
        for (i=1; i<n_points; i++) {
                Re_u_future[i] = Re_u_present[i]-2*(alpha*(Im_u_present[i+1]+Im_u_present[i-1])-2*(alpha+V[i]*dt)*Im_u_present[i]);
                Im_u_future[i] = Im_u_present[i]+2*(alpha*(Re_u_present[i+1]+Re_u_present[i-1])-2*(alpha+V[i]*dt)*Re_u_present[i]);
        
        //Variables to hold the previous value
        copy(Re_u_initial, Re_u_past, n_points);
        copy(Im_u_initial, Im_u_past, n_points);
        //Variables to hold the present value
        copy(Re_u_future, Re_u_present, n_points);
        copy(Im_u_future, Im_u_present, n_points);
        
            if(j%500 == 0){
        rho[i] = pow(Re_u_future[i],2) + Im_u_future[i]*Im_u_past[i];
        fprintf(in2,"%f ", rho[i]);
            }
    }
        
        fprintf(in2,"\n");
    }
    
    fclose(in2);
}
