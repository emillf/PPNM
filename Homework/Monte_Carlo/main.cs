using static System.Math;
using static System.Console;
using System.IO;
using System;
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
		WriteLine($"The plot shows that the error does indeed scale approximately as 1/sqrt(N)\n");
                using (StreamWriter writer = new StreamWriter("Vals1.dat")){
                	writer.WriteLine("# npoints err reerr errscal piapx pire ");
			int Nmax = 25000;
                	for(int N=0;N<Nmax;N+=500){
				var inte = plainmc(unitcircle,a1,b1,N);
                        	double npoints = N ;
                        	double err = inte.Item2;
				double errscal = 1/Sqrt(N);
				double reerr = Abs(err-PI);
				double piapx = inte.Item1;
				double pire = PI;
                        	writer.WriteLine($"{npoints} {err} {reerr} {errscal} {piapx} {pire}");
                        	}
                	}
		return 0;
	}//Main
}//Program
