using System;
using System.IO;
using static System.Math;
using static System.Console;
public class Program{
                public static matrix analytic_lagrange_jacobian(matrix A,vector x){
		//method to find analytic jacobian of f(v,λ)={Av-λv,v^Tv-1} for a given matrix A, vector v, and double λ where x={v,λ}
			int n = x.size-2;
			double λ = x[x.size-1];
                        matrix J=new matrix(n+2,n+2);
                        for(int j=0;j <= n;j++){ //Creating jacobian entries for i,j<n
                                for(int i=0; i<j+1 ;i++){
			 		if(i==n) J[i,j]=A[i,j]-λ; // we are in i,j<n case always because of loop conditions
					else J[i,j]=A[i,j];
					J[j,i]=J[i,j]; //Jacobian is symmetric in the i,j <= n case because A is symmetric
					}
                                }
			for(int i=0;i<=n;i++){ // We set λ deriviative to always appear in entries J[i,n+1]
				J[i,n+1]=-x[i];
				}
			for(int j=0;j<=n;j++){ //Handling f_n+1 entries without λ derivative
				J[n+1,j]=2*x[j];
				}
			J[n+1,n+1]=(double)0; // df_n+1/dλ case
                        return J;
                        }
		public static Func<vector,vector> Create_func_from_A(matrix A){ //Method to create Func<vector,vector> from A
			return (vector x) =>{
				vector v = new vector(x.size-1);
				for(int i=0;i<v.size;i++){
					v[i]=x[i];
					}
				double λ = x[x.size-1]; // My convention is still to take the last entry of x to be λ
				vector f_i = A*v-λ*v;
				vector fvals= new vector(x.size);
				for(int i=0; i<v.size; i++){
					fvals[i]=f_i[i];
					}
				fvals[x.size-1]= v.dot(v)-1;
				return fvals;
				};
			}
                public static (vector,double) newton_eigenvaluefinder(
                vector startv       // the start point for v
		,double startλ       // Startvalue for our lagrangian multiplier λ
		,matrix A            // Symmetric matrix A for to get the analytic Jacobian and our function to find roots of
                ,double acc=1e-3     // accuracy goal: on exit ‖f(x)‖ should be <acc
                ){
                        vector xnoλ=startv.copy();
			vector x = new vector(startv.size+1);
			for(int i=0;i<xnoλ.size;i++) x[i]=xnoλ[i];
			x[startv.size]=startλ;
			Func<vector,vector> f = Create_func_from_A(A);
                        vector fx=f(x),z,fz;
                        do{ /* Newton's iterations */
                                if(fx.norm() < acc) break; /* job done */
                                matrix J=analytic_lagrange_jacobian(A,x);
                                var QRJ = matrix.QR.decomp(J);
                                vector Dx = matrix.QR.solve(QRJ.Item1,QRJ.Item2,-fx); /* Newton's step */
                                double scaler=1.0;
                                double scalermin =Pow(10,-2);
                                do{ /* linesearch */
                                        z=x+scaler*Dx;
                                        fz=f(z);
                                        if( fz.norm() < (1-scaler/2)*fx.norm() ) break;
                                        if( scaler < scalermin ) break;
                                        scaler/=2;
                                        }while(true);
                                x=z; fx=fz;
                                }while(true);
			vector v= new vector(x.size-1);
			for(int i = 0; i<v.size; i++) v[i]=x[i];
			double λ = x[x.size-1];
                        return (v,λ);
                        }
	static int Main(){
		WriteLine("Part A) \n \nWe generate a random 3x3 symmetric matrix A\n");
		var rnd = new Random();
		matrix A = new matrix(3,3);
		for(int i=0; i<3; i++){
			for(int j=0; j<i+1;j++){
				double entry_ij=rnd.NextDouble()*10.0;
				A[i,j]=entry_ij;
				A[j,i]=entry_ij;
				}
			}
		A.print("A = \n");
		vector vstart= new vector(1.0,2.0,3.0);
		double λstart = 2.0;
		(vector v,double λ) = newton_eigenvaluefinder(vstart, λstart, A);
		WriteLine("\nUsing our method with initial guesses v=(1,2,3) and λ=2 we get:\n");
		v.print("v =");
		WriteLine($"\nλ =  {λ}");
		WriteLine("\nIf we calculate Av and λ*v we get\n");
		var λtest1= A*v;
		var λtest2= λ*v;
		λtest1.print("Av = ");
		λtest2.print("\nλv =");
		WriteLine("\nSo the method works as Av = λv");
		return 0;
		}//Main
	}//Program
