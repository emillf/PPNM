public class Vector{
   private double[] data;
   public int size => data.Length;
   public double this[int i]{   // indexer
	get => data[i];        // getter
	set => data[i]=value;  // setter
      }
	public Vector(int n){        // constructor
	data=new double[n];
      }
	public double Norm(Vector){
		double sumvec = 0;
		for(i=0;i<.size;i++){
			sumvec+=Vector[i]**2;
			}
		return System.Math.Sqrt(sumvec);
		}
	}
public class Matrix{
	public readonly int size1,size2;
	private double[] data;  // to keep matrix elements
	public matrix(int n,int m){      // constructor
		size1=n; size2=m;
		data = new double[size1*size2];
		}
	public double this[int i,int j]{     // indexer
		get => data[i+j*size1];
		set => data[i+j*size1]=value;
		}
	public matrix Transpose(matrix B){
		matrix T = B.copy()
		for(i=0;i<B.size2;i++){
			for(j=0;j<B.size1;j++){
				T[i,j]=B[j,i];
			}
		return T;
		}
	public static class QR{
		public static (matrix,matrix) decomp(matrix A){
			matrix Q=A.copy();
      			matrix R=new matrix(A.size2,A.size2);
			for(int i=0;i<A.size2; i++){
				R[i,i]=Q[i].norm();
				Q[i]/=R[i,i];
				for(int j =i+1;j<A.size2;j++){
					R[i,j]=Q[i].dot(Q(j));
					Q[j]-=Q[i]*R[i,j];
					}
				}
			return (Q,R);
      			}
   		public static Vector solve(matrix Q, matrix R, vector b){
				matrix U = Q.Transpose
				for (int i=b.size −1; i >=0; i−−){
				double sum=0;
					for (int k=i +1; k<b.size; k++){
						sum+=U[i,k]*b[k];
						b[i]=(b[i]−sum)/U[i,i];
					}
				}
			return b;
			}
   		public static double det(matrix R){
			double diag=0;
			double antidiag=0;
			for(int i=0; i<A.size2;i++){
				diag*=R[i,i];
				antidiag*=R[A.size2-i,A.size2-i];
				}
			return diag-antidiag;
			}

   		//public static matrix inverse(matrix Q,matrix R){
		//	}
		}
	}
