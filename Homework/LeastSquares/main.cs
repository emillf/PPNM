using static System.Math;
using System;
using static System.Console;
using static matrix;
using System.IO;
class Program{
	static (vector,matrix) lsfit(Func<double,double>[] fs, vector x, vector y, vector dy){
		int n = x.size;
		int m = fs.Length;
		var A = new matrix(n,m);
		var b = new vector(n);
		for(int i=0;i<n;i++){
			b[i]=y[i]/dy[i];
			for(int k=0;k<m;k++){
				A[i,k] = fs[k](x[i])/dy[i];
				}
			}
		var Asols = matrix.QR.decomp(A);
		var AQ=Asols.Item1;
		var AR=Asols.Item2;
		vector c = matrix.QR.solve(AQ,AR,b);
		var stupidsols = matrix.QR.decomp(A.T*A); // This is unnecessary and inefficient but it works and i dont have to change anything in my matrix class
		var ATAQ = stupidsols.Item1;
		var ATAR = stupidsols.Item2;
		var ATAI = matrix.QR.inverse(ATAQ,ATAR);
		matrix Sig = ATAI;
		return (c,Sig);
		}
	static int Main(){
		var rnd = new Random(4);
                var A = matrix.random(7,3,rnd);
                A.print("Random tall matrix A");
                WriteLine($" A decomposed into Q and R gives: \n");
                var QRtup =QR.decomp(A);
                matrix Q = QRtup.Item1;
                matrix R = QRtup.Item2;
                Q.print("Q =");
                R.print("R =");
                bool AeqQR=A.approx(Q*R);
                matrix prod = Q*R;
                prod.print("QR =");
                WriteLine($"QR = A? {AeqQR} \n");
		WriteLine("We see that R is upper triangular also, so it works as it should");
		var fs = new Func<double,double>[] { z => 1.0, z => -z};
		double[] ts = new double[] { 1,  2,  3, 4, 6, 9,   10,  13,  15};
		double[] acts = new double[] {117,100,88,72,53,29.5,25.2,15.2,11.1};
		double[] dys = new double[] {6,5,4,4,4,3,3,2,2};
		double[] lnys = new double[dys.Length];
		double[] lndys = new double[dys.Length];
		for(int i=0;i<ts.Length;i++){
			lnys[i] = Log(acts[i]);
			lndys[i] = dys[i]/acts[i];
			}
		vector t = new vector(ts);
                vector y = new vector(lnys);
                vector dy = new vector(lndys);
		var fit = lsfit(fs,t,y,dy);
		var c = fit.Item1;
		var a = Exp(c[0]);
		var lambda = c[1];
		var T_12 = Log(2)/lambda;
		WriteLine($"We get a = {a} and lambda={lambda}\n");
		WriteLine($"Half-life is T_1/2 = Log(2)/lambda = {T_12} days\n");
		WriteLine("Modern value is T_1/2 = 3.631(2) days\n");
		var S = fit.Item2;
		S.print("The covariance matrix:\n");
		var lambdaErr = Sqrt(S[1,1]);
		var aerr = Sqrt(S[0,0])/a;
		var T_12Err = lambdaErr * Log(2)/(Pow(lambda,2));
		WriteLine($"With use of covariance matrix and error propagation, the error on T_1/2 is found as {T_12Err}");
		WriteLine($"\nThis gives T_1/2 = {T_12} +- {T_12Err} days which is a bit outside the modern value {aerr}");
		int npoints = 100;
		double interval = 15-1;
		var Actfunc = new Func<double,double>(z=>a*Exp(-lambda*z));
		var Actfunclambda = new Func<double,double,double,double>((lus,aus,z)=>(a+aus)*Exp(-(lambda+lus)*z));
		using (StreamWriter writer = new StreamWriter("Vals.dat")){
                writer.WriteLine("# tss Activity texp yexp yexperr");
			for(int i=0;i<ts.Length;i++){
				double texp = ts[i];
				double yexp = acts[i];
				double yexperr = dys[i];
				writer.WriteLine($"{double.NaN} {double.NaN} {texp} {yexp} {yexperr}");
				}
                        for(int j=0;j<npoints;j++){
				double tss = interval/npoints*j;
                                double Activity=Actfunc(tss);
                                writer.WriteLine($"{tss} {Activity} {double.NaN} {double.NaN} {double.NaN}");
                                }
                        }
 		using (StreamWriter writer = new StreamWriter("OtherVals.dat")){
                writer.WriteLine("# tss1 lambdaplusaplus lambdaminusaminus lambdaplusaminus lambdaminusaplus");
                        for(int j=0;j<npoints;j++){
                                double tss1 = interval/npoints*j;
                                double lambdaplusaplus=Actfunclambda(lambdaErr,aerr,tss1);
				double lambdaminusaminus=Actfunclambda(-lambdaErr,-aerr,tss1);
				double lambdaplusaminus =Actfunclambda(lambdaErr,-aerr,tss1);
				double lambdaminusaplus = Actfunclambda(-lambdaErr,aerr,tss1);
                                writer.WriteLine($"{tss1} {lambdaplusaplus} {lambdaminusaminus} {lambdaplusaminus} {lambdaminusaplus}");
                                }
                        }
		return 0;
	}//Main
}//Program
