/* Project: 905 Project 1
 * File:   main.cpp
 * Author: alaina ross
 *
 * Created on January 12, 2017, 4:16 PM
 */

#include <iostream>
#include "../armadillo" //use this when using netbeans
//#include <armadillo> // need to change to this before turning in
using namespace std;
using namespace arma;

//prints vectors and matrix
void print_vals(mat , vec , vec ,int);

// performs gauss elimination
int main(int argc, char** argv) {
    int n;
    cout<<"n: ";
    cin >> n;
    
    mat a(n,n);
    vec b(n);
    vec x = zeros<vec>(n);
    
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        b(i)=i;
        for (int j=0;j<n;j++){
            a(i,j)=i+j+2;
            if(i==0 && j ==0)
                a(i,j)=1;
        }
    }
    
    cout<<"Before Gauss elimination"<<endl;
    print_vals(a,b,x,n);
    
    //do forward elimination
    double temp=0;
    for (int k=0;k<=n-2;k++){
        for (int i=k+1;i<=n-1;i++){
            temp=a(i,k)/a(k,k);
            for ( int j=k;j<=n-1;j++){
                a(i,j)=a(i,j)-a(k,j)*temp;
            }
            b(i)=b(i)-temp*b(k);
            if(abs(b(i))<=pow(10,-14))
                b(i)=0;
        }
    }

    //do back substitution
    if(abs(a(n-1,n-1))<=pow(10,-14)){
        x(n-1)=1;
        a(n-1,n-1)=0;
    }
    else{
        x(n-1)=b(n-1)/a(n-1,n-1);
    }
    double sum=0;
    for (int i=n-2;i>=0;i--){
        for (int j=n-1;j>i;j--){
            sum+=a(i,j)*x(j);
        }
        x(i)=(b(i)-sum)/a(i,i);
        if(abs(x(i))<=pow(10,-14))
            x(i)=0;
        sum=0;
    }
    
    cout<<"After Gauss elimination"<<endl;
    print_vals(a,b,x,n); 
    
    
    return 0;
}

void print_vals(mat A, vec b, vec x,int n){
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
    cout<<"x: ";
    for (int i=0;i<n;i++){
        cout<<x(i)<<" ";
    }
    cout<<endl;
    
}