/* Project: 905 Project 1
 * File:   main_general.cpp
 * Author: alaina ross
 *
 * Created on January 12, 2017, 4:16 PM
 */

#include <iostream>
#include "../armadillo" //use this when using netbeans
//#include <armadillo> // need to change to this before turning in
#include <sys/time.h>
using namespace std;
using namespace arma;

//template functions
void print_vals(mat , vec , vec ,int);
void get_walltime_( double* wcTime);
void get_walltime( double* wcTime);

// performs gauss elimination
// for general matrix
int main(int argc, char** argv) {
    int n;
    cout<<"n: ";
    cin >> n;
    double h = 1/(double(n)+1);
    double start=0;
    double end=0;
    
    mat a = zeros<mat>(n,n);
    vec b(n);
    vec v = zeros<vec>(n);
    vec x(n);
    
    //initialize x values
    x(0)=h;
    for (int i=1; i<n ;i++)
        x(i)=x(i-1)+h;

    
    //initialize matrix and vector
//    for (int i=0;i<n;i++){
//        b(i)=i;
//        for (int j=0;j<n;j++){
//            a(i,j)=i+j+2;
//            if(i==0 && j ==0)
//                a(i,j)=1;
//        }
//    }
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j)
                a(i,j)=2;
            else if (i==j+1 or i==j-1)
                a(i,j)=-1;
        }
        b(i)=h*h*100*exp(-10*(x(i)));
    }
    
    if(n<=10){
        cout<<"Before Gauss elimination"<<endl;
        print_vals(a,b,v,n);
    }
    
    get_walltime(&start);
    
    //do forward elimination
    double temp=0;
    for (int k=0;k<=n-2;k++){
        for (int i=k+1;i<=n-1;i++){
            temp=a(i,k)/a(k,k);
            for ( int j=k;j<=n-1;j++){
                a(i,j)-=a(k,j)*temp;
            }
            b(i)-=temp*b(k);
            if(abs(b(i))<=pow(10,-14))
                b(i)=0;
        }
    }

    //do back substitution
    if(abs(a(n-1,n-1))<=pow(10,-14)){
        v(n-1)=1;
        a(n-1,n-1)=0;
    }
    else{
        v(n-1)=b(n-1)/a(n-1,n-1);
    }
    double sum=0;
    for (int i=n-2;i>=0;i--){
        for (int j=n-1;j>i;j--){
            sum+=a(i,j)*v(j);
        }
        v(i)=(b(i)-sum)/a(i,i);
        if(abs(v(i))<=pow(10,-14))
            v(i)=0;
        sum=0;
    }
    
    get_walltime(&end);
    
    if(n<=10){
        cout<<"After Gauss elimination"<<endl;
        print_vals(a,b,v,n); 
    }
    
    cout<<"Total walltime (sec) : "<<end-start<<endl;
    return 0;
}

void print_vals(mat A, vec b, vec v,int n){
    cout<<"A: ";
    for (int i=0;i<n;i++){
        if(i>0){
            cout<<"   ";
         }
        for (int j=0;j<n;j++){
            cout<<A(i,j)<<" ";
        }
        cout<<endl;
    }
    cout<<"b: ";
    for (int i=0;i<n;i++){
        cout<<b(i)<<" ";
    }
    cout<<endl;
    cout<<"v: ";
    for (int i=0;i<n;i++){
        cout<<v(i)<<" ";
    }
    cout<<endl;
    
}

void get_walltime_( double* wcTime) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  *wcTime = ( double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void get_walltime( double* wcTime) {
  get_walltime_(wcTime);
}