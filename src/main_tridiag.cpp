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
void print_vals(vec,vec, vec, vec, vec,int);

// performs gauss elimination
int main(int argc, char** argv) {
    int n;
    cout<<"n: ";
    cin >> n;
    double h = 1/(double(n)+1);
    
    vec b(n);
    vec v = zeros<vec>(n);
    vec a(n-1);
    vec c(n-1);
    vec d(n);
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
        d(i)=2;
        b(i)=h*h*100*exp(-10*(x(i)));
        if (i<n-1){
            a(i)=-1;
            c(i)=-1;
        }
    }
    
    if(n<=10){
        cout<<"Before Gauss elimination"<<endl;
        print_vals(a,c,d,b,v,n);            
    }
    
    //do forward elimination
    double temp=0;
    for (int k=0;k<=n-2;k++){
        temp=a(k)/d(k);
        a(k)=0;
        d(k+1)-=c(k)*temp;
        b(k+1)-=temp*b(k);
    }

    //do back substitution
    if(abs(d(n-1))<=pow(10,-14)){
        v(n-1)=1;
        d(n-1)=0;
    }
    else{
        v(n-1)=b(n-1)/d(n-1);
    }
    for (int i=n-2;i>=0;i--){
        v(i)=(b(i)-c(i)*v(i+1))/d(i);
    }

    
    if(n<=10){
        cout<<"After Gauss elimination"<<endl;
        print_vals(a,c,d,b,v,n);
    }
    
    
    return 0;
}

void print_vals(vec a,vec c, vec d, vec b, vec v,int n){
    cout<<"A: "<<d(0)<<" "<<c(0)<<" ";
    for (int j=3;j<=n;j++)
        cout<<"0 ";
    cout<<endl;
    for (int i=1;i<n;i++){
        cout<<"   ";
        for (int j=0;j<i-1;j++)
            cout<<"0 ";
        cout<<a(i-1)<<" "<<d(i)<<" ";
        if(i<n-1)
            cout<<c(i)<<" ";
        for (int j=i+2;j<n;j++)
            cout<<"0 ";
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