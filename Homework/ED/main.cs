using System;
using static System.Math;
public static class jacobi{
	public static void timesJ(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int i=0;i<A.size1;i++){
			double aip = A[i,p];
			double aiq = A[i,q];
			A[i,p]=c*aip-s*aiq;
			A[i,q]=s*aip+c*aiq;
			}
		}
	public static void Jtimes(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int j =0;j<A.size2;j++){
			double apj = A[p,j];
			double aqj = A[q,j];
			A[p,j]=c*apj-s*aqj;
			A[q,j]=s*apj+c*aqj;
			}
		}

	public static (vector,matrix) cyclic(matrix M){
		matrix A=M.copy();
		matrix V=matrix.id(M.size1);
		vector w=new vector(M.size1);
		for(int i =
		/* run Jacobi rotations on A and update V */
		/* copy diagonal elements into w */
		return (w,V);
		}
	}
