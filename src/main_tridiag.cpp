/* Project: 905 Project 1
 * File:   main_tridiag.cpp
 * Author: alaina ross
 *
 * Created on January 24, 2017, 1:45 PM
 */

#include <iostream>
#include<fstream>
#include <armadillo>
#include<time.h>
using namespace std;
using namespace arma;

//prints vectors and matrix
void print_vals(vec,vec, vec, vec, vec,int);

// performs gauss elimination
// for general tridiagonal matrix
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
    clock_t begin,finish;   
    
    //initialize x values
    x(0)=h;
    for (int i=1; i<n ;i++)
        x(i)=x(i-1)+h;

    //initialize vectors
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
    
    begin=clock();    
    
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
    
    finish=clock();

    
    if(n<=10){
        cout<<"After Gauss elimination"<<endl;
        print_vals(a,c,d,b,v,n);
    }
    
    cout<<"Total CPU time (sec) : "<<((double)finish-(double)begin)/
            CLOCKS_PER_SEC<<endl;
    
    vec u(n);
    ofstream outfile("output_tridiag.txt");
    ofstream error("error_tridiag.txt");
    ofstream actual("output_analytic.txt");
    outfile<<"0.0 0.0"<<endl;
    actual<<"0.0 0.0"<<endl;
    error<<log10(h)<<endl;
    for (int i=0;i<n;i++){
        outfile<<x(i)<<" "<<v(i)<<endl;
        u(i)=1-(1-exp(-10))*x(i)-exp(-10*x(i));
        actual<<x(i)<<" "<<u(i)<<endl;
        error<<x(i)<<" "<<log10(abs((v(i)-u(i))/u(i)))<<endl;
    }
    outfile<<"1.0 0.0"<<endl;
    actual<<"1.0 0.0"<<endl;
    outfile.close();
    actual.close();
    error.close();
    
    
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