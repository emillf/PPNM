using System;
using System.IO;
using static System.Math;
using static System.Console;
using System.Diagnostics;
public class Program{
                public static matrix analytic_lagrange_jacobian(matrix A,vector x){
		//method to find analytic jacobian of f(v,λ)={Av-λv,v^Tv-1} for a given matrix A, vector v, and double λ where x={v,λ}
			int n = x.size-2;
			double λ = x[x.size-1];
                        matrix J=new matrix(n+2,n+2);
                        for(int j=0;j <= n;j++){ //Creating jacobian entries for i,j<n
                                for(int i=0; i<j+1 ;i++){
			 		if(i==j) J[i,j]=A[i,j]-λ; // we are in i,j<n case always because of loop conditions
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
			J[n+1,n+1]=0.0; // df_n+1/dλ case
                        return J;
                        }
		public static Func<vector,vector> Create_func_from_A(matrix A){ //Method to create Func<vector,vector> from symmetric A
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
		public static matrix Create_H(double delta_r =0.05,double rmin=0.1, double rmax=20.0,double l=0.0){
			int n = (int)((rmax - rmin)/delta_r);
			matrix H = new matrix(n,n);
			double coeff= 1.0/(delta_r*delta_r);
			for(int i=0;i<n;i++){ //We drop H[n,n] since we set u(rmin)=u(rmax)=0
				double ri =rmin + (i)*delta_r;
				H[i,i]=-coeff-1.0/ri + l*(l+1)/(2*ri*ri);
				if(i>0)H[i,i-1]=0.5*coeff;
				if(i<n-1) H[i,i+1]=0.5*coeff;
				}
			//H[0, 0] = 0.0;
			//H[1,0] =0.0;
			//H[n - 1, n - 1] = 0.0;
			//H[n-2,n-1]=0.0; // We set first and last column to be 0 to enforce boundary conditions
			return H;
			}
                public static (vector,double,double) newton_eigenvaluefinder(
                vector startv       // the start point for v
		,double startλ       // Startvalue for our lagrangian multiplier λ
		,matrix A            // Symmetric matrix A for to get the analytic Jacobian and our function to find roots of
                ,double acc=1e-3     // accuracy goal: on exit ‖f(x)‖ should be <acc
		,bool analytic=false //Option to use analytic linesearch instead of backtracking linesearch
                ){
			var stopwatch=new Stopwatch();
			int maxiterations =(int)Pow(10,5);
			stopwatch.Start();
                        vector xnoλ=startv.copy();
			vector x = new vector(startv.size+1);
			for(int i=0;i<xnoλ.size;i++) x[i]=xnoλ[i];
			x[startv.size]=startλ;
			Func<vector,vector> f = Create_func_from_A(A);
                        vector fx=f(x),z,fz;
			int iteration=0;
                        do{ /* Newton's iterations */
                                if(fx.norm() < acc) break; /* job done */
        			if (iteration >= maxiterations) {
            				WriteLine($"Newton method reached max iterations ({maxiterations}) without converging.");
            				break;
					}
				if (double.IsNaN(fx.norm()) || fx.norm() > 1e10) {
					WriteLine("Divergence detected.");
    					break;
					}
                                matrix J=analytic_lagrange_jacobian(A,x);
                                var QRJ = matrix.QR.decomp(J);
                                vector Dx = matrix.QR.solve(QRJ.Item1,QRJ.Item2,-fx); /* Newton's step */
                                double scaler=1.0;
                                double scalermin =Pow(2,-6);
				double alpha= 1.0; //Placeholder for optimal analytic scaler
				if(analytic){ //Analytic linesearch
					vector Jdx=J*Dx;
					alpha=-fx.dot(Jdx)/(Jdx.dot(Jdx));
					z=x+alpha*Dx;
					fz=f(z);
					}
				else{
                                	do{ // Numerical linesearch
                                        	z=x+scaler*Dx;
                                        	fz=f(z);
                                        	if( fz.norm() < (1-scaler/2)*fx.norm() ) break;
                                        	if( scaler < scalermin ) break;
                                        	scaler/=2;
                                        	}while(true);
					}
                                	x=z; fx=fz; iteration++;
                                }while(true);
			stopwatch.Stop();
			vector v= new vector(x.size-1);
			for(int i = 0; i<v.size; i++) v[i]=x[i];
			double λ = x[x.size-1];
                        return (v,λ,stopwatch.Elapsed.TotalSeconds);
                        }
	static int Main(){
		WriteLine("Part A) \n \nWe generate a random 3x3 symmetric matrix A with values from 0 to 10\n");
		var rnd = new Random();
		int ndim=6;
		matrix A = new matrix(ndim,ndim);
		for(int i=0; i<ndim; i++){
			for(int j=0; j<i+1;j++){
				double entry_ij=rnd.NextDouble()*10.0;
				A[i,j]=entry_ij;
				A[j,i]=entry_ij;
				}
			}
		A.print("A = \n");
		vector vstart= new vector(ndim);
		vector vstarttest = new vector(ndim+1);
		for(int i =0;i<ndim;i++){
			vstart[i]=1.0;
			vstarttest[i]=1.0;
			}
		vstarttest[ndim]=2.0;
		double λstart = 2.0;
		var Jacobian_test = analytic_lagrange_jacobian(A,vstarttest);
		Jacobian_test.print("\n");
		vstarttest.print("\n");
		(vector v,double λ,double time1) = newton_eigenvaluefinder(vstart, λstart, A);
		WriteLine("\nUsing our method with initial guesses v = (1,1,1) and λ = 2 we get:\n");
		v.print("v =");
		WriteLine($"\nλ =  {λ}");
		WriteLine($"\nIf we calculate Av and λ*v we get in time {time1} s\n");
		var λtest1= A*v;
		var λtest2= λ*v;
		λtest1.print("Av = ");
		λtest2.print("\nλv =");
		WriteLine("\nSo the method works as we see Av = λv");
		WriteLine("\nWe time the process for a larger and larger n x n symmetric matrix in the plot.");
		WriteLine("\nTo remove variation in the quality of the starting guess, the n x n matrix will be a tridiagonal toeplitz matrix.");
		WriteLine("\nEntries are defined as A[i,i] = -1, and A[i+1,i]=A[i,i+1]=2 (for 1<=i<=n-1),");
		WriteLine("\nWith the kth eigenvalue λ_k = 2-2*cos(k*pi/(n+1))");
		WriteLine("\nand corresponding eigenvector entries v_j=sqrt(2/(n+1))*sin(j*k*pi/(n+1))");
		WriteLine("\nThis ensures well spaced eigenvalues and some offdiagonal terms. \nI can thus choose λstart = 2 and vstart=(1,1,1,1,....)");
		WriteLine("\nAs can be seen in Time_plot time goes like O(n^2), so generally around n=400 it becomes unfeasible for my box");
		WriteLine("\nThe exact time it becomes unfeasible also depends on the matrix choice and guess quality of course.");
		WriteLine($"\nI tried a random 100x100 matrix with the same start guesses as the toeplitz this took about twice as long as the 100x100 toeplitz");
		int nmax=20;
		int index=0;
		vector timetoes = new vector(10);
		for(int n = nmax/10;n<nmax+nmax/10;n+=nmax/10){
		matrix Atoe=new matrix(n,n);
		for(int i=0;i<n;i++){
			Atoe[i,i]=2.0;
			}
		for(int j = 1; j<n;j++){
			Atoe[j-1,j]=-1.0;
			Atoe[j,j-1]=-1.0;
			}
                for(int j = 0; j<n-1;j++){
                        Atoe[j+1,j]=-1.0;
                        Atoe[j,j+1]=-1.0;
                        }
		double λstarttoe=2.0;
		vector vstarttoe= new vector(Atoe.size1);
		vector vtesttoe= new vector(Atoe.size1);
		for(int i=0; i<Atoe.size1; i++){
			vstarttoe[i] = 1;
			vtesttoe[i] = Sqrt(2/(Atoe.size1))*Sin(((i+1)*1*PI)/(n+1));
			}
		(vector vtoe, double λtoe, double timetoe)=newton_eigenvaluefinder(vstarttoe,λstarttoe,Atoe);
		timetoes[index]=timetoe;
		index++;
		}
		double rmax=30.0;
		double rmin =0.05;
		double delta_r = 0.05;
		matrix Hl0=Create_H(delta_r,rmin,rmax);
		matrix Hl1=Create_H(delta_r,rmin,rmax);
		int H_n = Hl0.size1;
		vector ustartn1 = new vector(H_n);
		vector ustartn2l0= new vector(H_n);
		vector ustartn2l1= new vector(H_n);
		for (int i = 0; i < H_n; i++) {
    			double r =rmin+ i*delta_r;
    			ustartn1[i] = 2.0*r*Exp(-r)+0.1;
			}
                for (int i = 0; i < H_n; i++) {
                	double r =rmin+ i*delta_r;
                        ustartn2l0[i] = (r/Sqrt(2))*(1.0-r/2.0)*Exp(-r/2.0)+0.1;
                         }
                for (int i = 0; i < H_n; i++) {
                	double r =rmin+ i *delta_r;
                        ustartn2l1[i] = 1/Sqrt(24)*r*r * Exp(-r/2.0)+0.1;
                        }
		(vector un1,double E1,double timen1) = newton_eigenvaluefinder(ustartn1,-1/2.1,Hl0,1e-3,true);
		(vector un2l0,double E2l0,double timen2l0) = newton_eigenvaluefinder(ustartn2l0,-1/8.1,Hl0,1e-3,true);
		(vector un2l1,double E2l1,double timen2l1) =  newton_eigenvaluefinder(ustartn2l0,-1/8.1,Hl1,1e-3,true);
		double utestn1 = (Hl0*un1).dot(Hl0*un1)-(E1*un1).dot(E1*un1);
		double utestn2l0 = (Hl0*un2l0).dot(Hl0*un2l0)-(E2l0*un2l0).dot(E2l0*un2l0);
		double utestn2l1 = (Hl1*un2l1).dot(Hl1*un2l1)-(E2l1*un2l1).dot(E2l1*un2l1);
		WriteLine($"\n\nPart C)\n\n Calculating several lowest eigenfunctions of the hydrogen atom requires solving the eigenvalue problem Hu=Eu");
		WriteLine("\nWhere H[i,i]=-1/h^2+l(l+1)/(2r[i]^2)-1/r_i and H[i,i±1]=1/(2h^2), where h is small.");
		WriteLine("\nThis equation has been derived by expanding the derivatives of u the radial Schrödinger equation in atomic units.");
		WriteLine($"\nLike in our roots exercise, we choose rmax to be {rmax} and rmin={rmin}");
		WriteLine("\nWe calculate the n=1 l=0, n=2 l=0 and n=2 l=1  lowest states of the Hydrogen atom");
		WriteLine($"\n These give us the energies {E1} {E2l0} and {E2l1} respectively (in hartrees)");
		WriteLine($"\nAlso see that ||H_(l=0)*u_10||^2 - ||E_1*un1||^2={utestn1:F3},||H_(l=0)*u_20||^2 - ||E_2*u_20||^2={utestn2l0:F3}");
		WriteLine($"\nand that |H_(l=1)*u_21||^2 - ||E2*u_21||^2={utestn2l0:F3}");
		WriteLine("\nSo the wrong energies (they should be -1/2,-1/8,-1/8) are not due to the linesearch failing.");
		WriteLine("\nRather it is likely due to needing better step length, smaller rmin and rmax which is not feasible as H becomes too large,");
		WriteLine("\nOr a way to write H which enforces the boundary conditions better - namely u(rmin)=rmin-rmin^2");
		WriteLine("\nA plot of the corresponding reduced radial wavefunctions can be found in Eigenfuncs.svg");
                using (StreamWriter writer = new StreamWriter("Times.dat")){
                	writer.WriteLine("#ns ts");
                	for(int i=1;i<10;i++){
                                var ns = (double)nmax*i/10;
                                var ts = timetoes[i-1];
                        	writer.WriteLine($"{ns} {ts}");
                                }
                        }
		using (StreamWriter writer = new StreamWriter("Eigenfuncs.dat")){
			writer.WriteLine("#r #un1 #un2l0 #un2l1 #urealn1 #urealn2l0 #urealn2l1");
			for(int i =0; i< un1.size; i++){
				var r = rmin + (i)*delta_r;
				var un1s=un1[i];
				var un2l0s= un2l0[i];
				var un2l1s= un2l1[i];
				var un1r = 2.0*r*Exp(-r);
				var un2l0r = (r/Sqrt(2))*(1.0-r/2.0)*Exp(-r/2.0);
				var un2l1r = 1/Sqrt(24)*r*r * Exp(-r/2.0);
				writer.WriteLine($"{r} {un1s} {un2l0s} {un2l1s} {un1r} {un2l0r} {un2l1r}");
				}
			}
		return 0;
		}//Main
	}//Program
