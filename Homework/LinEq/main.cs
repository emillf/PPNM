using System;
using static System.Console;
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
                double antidiag=1;
                for(int i=0; i<R.size2;i++){
                        diag*=R[i,i];
                        antidiag*=R[R.size2-i,R.size2-i];
                        }
                return diag-antidiag;
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
public class Program{
	static void Main(){
        	var rnd = new Random(3);
		var A = matrix.random(7,3,rnd);
		A.print("Random matrix A");
		WriteLine($" A decomposed into Q and R gives: /n");
		var QRtup =QR.decomp(A);
		matrix Q = QRtup.Item1;
		matrix R = QRtup.Item2;
		Q.print("Q =");
		R.print("R =");
		bool AeqQR=A.approx(Q*R);
		matrix prod = Q*R;
		prod.print("QR =");
		WriteLine($"QR = A? {AeqQR} /n");
		matrix QTQ = Q.T*Q;
		QTQ.print("Q^T Q =");
		WriteLine("Now lets check if solve works");
		WriteLine("Our random symmetric matrix and vector respectively:");
		var B = matrix.random(5,5,rnd);
		var b = vector.random(5,rnd);
		B.print("B = ");
		b.print("b = ");
		WriteLine("Solving Ax=b using QR decomp gives us");
		var QRtupnew = QR.decomp(B);
		matrix Qnew = QRtupnew.Item1;
		matrix Rnew = QRtupnew.Item2;
		Qnew.print("Q = ");
		Rnew.print("R = ");
		vector x = QR.solve(Qnew,Rnew,b);
		WriteLine("Solving QRx=b we get: x =");
		x.print();
		var Bxprod =B*x;
		WriteLine("With Bx =");
		Bxprod.print();
		bool Bxeqb = b.approx(Bxprod);
		WriteLine($"Bx = b ? {Bxeqb}");
		}
	}
