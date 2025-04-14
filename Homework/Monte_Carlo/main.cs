using static System.Math;
using static System.Console;
using System.IO;
using System;
using System.Collections.Generic;
public class Program{
        static (double,double) plainmc(Func<vector,double> f,vector a,vector b,int N){
                int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
                double sum=0,sum2=0;
                var x=new vector(dim);
                var rnd=new Random();
                for(int i=0;i<N;i++){
                        for(int k=0;k<dim;k++)x[k]=a[k]+rnd.NextDouble()*(b[k]-a[k]);
                        double fx=f(x); sum+=fx; sum2+=fx*fx;
                                }
                double mean=sum/N, sigma=Sqrt(sum2/N-mean*mean);
                var result=(mean*V,sigma*V/Sqrt(N));
                return result;
                }
	static double corput(int n, int b){
		double q = 0;
		double bk = (double)1/b;
		while(n>0){
			q+=(n%b)*bk;
			n/=b;
			bk/=b;
			}
		return q;
		}
	static List<int> primelist(int n){
		List<int> primes = new List<int>();
		int candidate = 2;
		for(int i=0;i<n;){
			bool IsPrime=true;
			for(int j=0; j<primes.Count; j++){
				if(primes[j]*primes[j]>candidate) break;
				if(candidate%primes[j]==0){
					IsPrime= false;
					break;
					}
				} // nested for loop
			if(IsPrime==true){
				primes.Add(candidate);
				i++;
				}
			candidate++;
			}//parent for loop
		return primes;
		} //method
	static void halton(int n, vector x){
		List<int> bases = primelist(x.size);
		for(int i = 0; i<x.size; i++){
			x[i] = corput(n,bases[i]);
			}
                }
        static void lattice(int n, vector x){
                List<double> alphas = new List<double>();
                List<int> primes = primelist(x.size);
                for(int i = 0;i<x.size;i++)alphas.Add(Sqrt(primes[i])%1);
                for(int i = 0;i<x.size;i++)x[i]=(n*alphas[i]%1);
                }
	static (double,double) quasimc(Func<vector,double> f, vector a, vector b, int N, bool halt=true){
		int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
                double sum=0,sum2=0;
                var x = new vector(dim);
                var haltonvec=new vector(dim);
                var latticevec = new vector(dim);
                for(int i=0;i<N;i++){
			halton(i,haltonvec);
			for(int j = 0; j<dim; j++)x[j]=a[j]+(b[j]-a[j])*haltonvec[j];
                        sum+=f(x);
                        lattice(i,latticevec);
                        for(int k = 0; k<dim; k++)x[k]=a[k]+(b[k]-a[k])*latticevec[k];
                        sum2+=f(x);
                        }
		double meanhalt=sum/N;
                double meanlat=sum2/N;
                double sigma = Abs(meanhalt-meanlat);
                if(halt==true) return (meanhalt,sigma);
                else return (meanlat,sigma);
                }
	static int Main(){
		WriteLine($"Part a)\n We do integrals using our plain montecarlo integration:\n ");
		// We calculate the unit circle in polar coordinates
		Func<vector,double> unitcircle = v => v[0];
		vector a1 = new vector(0,0);
		vector b1 = new vector(1,2*PI);
		int npoints1 = 10000;
		var Aunitcircle = plainmc(unitcircle,a1,b1,npoints1);
                WriteLine($" First we check the area of the unit circle using {npoints1} points\n");
		WriteLine($"(from θ=0 to θ=2π and r=0 to r=1) ∫∫rdxdθ = {Aunitcircle.Item1}+-{Aunitcircle.Item2}\n");
		WriteLine($"The plot shows that the error scales very roughly as 1/sqrt(N)\n \n");
		int npoints2 = 1000000;
		WriteLine($"Now lets do the very difficult singular integral (here called I) using {npoints2} points\n");
                Func<vector,double> difffunc = v => 1/Pow(PI,3)*1/(1-Cos(v[0])*Cos(v[1])*Cos(v[2]));
		vector a2= new vector(0,0,0);
		vector b2= new vector(PI,PI,PI);
		var diffval = plainmc(difffunc,a2,b2,npoints2);
		WriteLine($"This gives I={diffval.Item1}+-{diffval.Item2}\n");
		WriteLine($"The real error is {1.39320392-diffval.Item1}\n");
		using (StreamWriter writer = new StreamWriter("Vals1.dat")){
                	writer.WriteLine("# npoints err reerr errscal piapx pire ");
			int Nmax = 25000;
                	for(int N=0;N<Nmax;N+=500){
				var inte = plainmc(unitcircle,a1,b1,N);
                        	double npoints = N ;
                        	double err = inte.Item2;
				double errscal = 1/Sqrt(N);
				double reerr = Abs(inte.Item1-PI);
				double piapx = inte.Item1;
				double pire = PI;
                        	writer.WriteLine($"{npoints} {err} {reerr} {errscal} {piapx} {pire}");
                        	}
                	}
		var diffvalq = quasimc(difffunc,a2,b2,npoints1);
		WriteLine($"Part B: \n as an example Calculating I with our (halton) quasirandom sampling using lattice as error estimate sequence using {npoints1} points we get I = {diffvalq.Item1}+-{diffvalq.Item2}");
		WriteLine("\n We check the error scaling of our quasirandom sampling on the same integral as our plain sampling in the plot \n");
                using (StreamWriter writer = new StreamWriter("Vals2.dat")){
                        writer.WriteLine("# npoints errplain rerrplain errscalplain errquasi reerrquasi");
                        int Nmax = 25000;
                        for(int N=0;N<Nmax;N+=500){
				var intequasi = quasimc(unitcircle,a1,b1,N);
				var inteplain = plainmc(unitcircle,a1,b1,N);
                                double npoints = N ;
				double errquasi = intequasi.Item2;
				double errplain = inteplain.Item2;
                                double errscalplain = 1/Sqrt(N);
                                double reerrquasi = Abs(intequasi.Item1-PI);
                                double reerrplain = Abs(inteplain.Item1-PI);
                                writer.WriteLine($"{npoints} {errplain} {reerrplain} {errscalplain} {errquasi} {reerrquasi}");
                                }
                        }
		return 0;
	}//Main
}//Program
