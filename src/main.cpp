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
        b(i)=i-1;
        for (int j=0;j<n;j++){
            a(i,j)=i+j+1;
            if(i==1 && j ==1)
                a(i,j)=1;
        }
    }
    
    double temp=0;
    for (int k=0;k<=n-2;k++){
        cout<<"1"<<endl;
        for (int i=k+1;i<=n-1;i++){
            cout<<"2"<<endl;
            temp=a(i,k)/a(k,k);
            for ( int j=k;j<=n-1;j++){
                cout<<"3"<<endl;
                a(i,j)=a(i,j)-a(k,j)*temp;
            }
            b(i)=b(i)-temp*b(k);
        }
    }
    
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