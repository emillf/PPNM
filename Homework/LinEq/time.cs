using System;
public static class QR{
	public static (matrix,matrix) decomp(matrix A){
                matrix Q=A.copy();
                matrix R=new matrix(A.size2,A.size2);
                for(int i=0;i<A.size2; i++){
                        R[i,i]=Q[i].norm();
                        Q[i]/=R[i,i];
                        for(int j =i+1;j<A.size2;j++){
                                R[i,j]=Q[i].dot(Q[j]);
                                Q[j]-=Q[i]*R[i,j];
                                }
                        }
                return (Q,R);
                }
        public static vector solve(matrix Q, matrix R, vector b){
                var c = Q.T*b;
                var U = R.copy();
                for(int i=c.size -1; i >=0; i--){
                        double sum=0;
                        for (int k=i +1; k<c.size; k++) sum+=U[i,k]*c[k];
                        c[i]=(c[i]-sum)/U[i,i];
                        }
                return c;
                }
        public static double det(matrix R){
                double diag=1;
                for(int i=0; i<R.size2;i++){
                        diag*=R[i,i];
                        }
 		return diag;
                }
        public static matrix inverse(matrix Q,matrix R){
           int n = Q.size1;
           matrix Ainv = new matrix(n,n);
           vector unit = new vector(n);
           for(int i=0;i<n;i++){
                   unit[i] = 1;
                   Ainv[i] = solve(Q,R,unit);
                   unit[i] = 0;
                }
           return Ainv;
        }
}
public class Time{
	static int Main(string[] args){
		var rnd = new Random(3);
                int N = 5;
		foreach(var arg in args){
			string[] words = arg.Split(":");
			if(words[0]=="-size") N = int.Parse(words[1]);
			}
		matrix A = matrix.random(N,N,rnd);
		var QRtup = QR.decomp(A);
		var Q = QRtup.Item1;
		var R = QRtup.Item2;
		Q.print();
		R.print();
		return 0;
		}
	}
