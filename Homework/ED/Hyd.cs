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

        public static (matrix M,vector w, matrix V) cyclic(matrix M){
                matrix A=M.copy();
                matrix V=matrix.id(M.size1);
                vector w=new vector(M.size1);
                bool changed;
                double tolerance = Pow(10,-18);
                int maxIterations = 100;
                int iteration = 0;
                do{
                        double[] diagOld = new double[A.size1];
                        for (int i = 0; i < A.size1; i++){
                                diagOld[i] = A[i, i];
                                }
                        changed = false;
                        for(int p = 0;p<A.size1-1;p++){
                                for(int q=p+1;q<A.size1;q++){
                                        double apq=A[p,q];
                                        double aqq=A[q,q];
                                        double app=A[p,p];
                                        if(Abs(apq)>tolerance)
                                                {
                                                double theta=0.5*Atan2(2*apq,aqq-app);
                                                timesJ(A,p,q,theta);
                                                Jtimes(A,p,q,theta);
                                                timesJ(V,p,q,theta);
                                                A[p,q]=0;
                                                A[q,p]=0;
                                                }
                                        }
                                }
                        for (int i = 0; i < A.size1; i++){
                                if (Abs(A[i, i]-diagOld[i])>tolerance)
                                        {
                                        changed = true;
                                        break;
                                        }
                                }
                iteration++;
                }while(changed && iteration<maxIterations);
                for(int i = 0 ;i<A.size1;i++){
                        w[i]=A[i,i];
                        }
                return (A,w,V);
                }
	}
class Program{
	static int Main(string[] args){
		rmax=10
		dr=0.3
		foreach(var arg in args){
			string[] words = arg.Split(":");
			if(words[0]=="-rmax") rmax = double.Parse(words[1]);
			if(word[2]=="-dr") dr = double.Parse(words[3];
                        }
		
