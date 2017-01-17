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
    
    mat A(n,n);
    vec b(n);
    vec x = zeros<vec>(n);
    
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        b(i)=i-1;
        for (int j=0;j<n;j++){
            A(i,j)=i+j+1;
            if(i==1 && j ==1)
                A(i,j)=1;
        }
    }
    
    print_vals(A,b,x,n);
    
    
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