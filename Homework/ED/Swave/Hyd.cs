using System;
using static System.Math;
using System.Collections.Generic;
using System.IO;
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
		double rmax=10;
		double dr=0.3;
		string mode = "both";
        	for (int i = 0; i < args.Length; i++){
            		if (args[i] == "-rmax" && i+1 < args.Length){
                		rmax = double.Parse(args[i+1]);
                		i++;
            			}
            		else if (args[i] == "-dr" && i+1 < args.Length){
                		dr = double.Parse(args[i+1]);
                		i++;
            			}
            	else if (args[i] == "-mode" && i+1 < args.Length){
                	mode = args[i+1].ToLower();
                	i++;
            		}
        	}
        Console.WriteLine($"Running with rmax={rmax}, dr={dr}, mode={mode}");
        if (mode == "dr" || mode == "both"){
            TestDrConvergence(rmax);
        }
        if (mode == "rmax" || mode == "both"){
            TestRmaxConvergence(dr);
        }
        return 0;
    }
    static void TestDrConvergence(double fixedRmax)
    {
        Console.WriteLine($"\nTesting dr convergence with fixed rmax = {fixedRmax}");
        // Create a list of dr values to test
        // Using a logarithmic distribution to cover different scales
        List<double> drValues = new List<double>();
        for (double dr = 0.1; dr <= 1.0; dr += 0.1)
        {
            drValues.Add(dr);
        }
        // Also add some finer values at the lower end
        drValues.Add(0.01);
        drValues.Add(0.02);
        drValues.Sort();
        // Calculate ground state energy for each dr
        List<double> energies = new List<double>();
        foreach (double dr in drValues)
        {
            double energy = CalculateGroundStateEnergy(fixedRmax, dr);
            energies.Add(energy);
            Console.WriteLine($"dr = {dr}, ε₀ = {energy}");
        }
        // Write results to a file for plotting
        using (StreamWriter writer = new StreamWriter("dr_convergence.dat"))
        {
            writer.WriteLine("# dr ε₀");
            for (int i = 0; i < drValues.Count; i++)
            {
                writer.WriteLine($"{drValues[i]} {energies[i]}");
            }
        }
        Console.WriteLine("Results saved to dr_convergence.dat");
    }
    static void TestRmaxConvergence(double fixedDr)
    {
        Console.WriteLine($"\nTesting rmax convergence with fixed dr = {fixedDr}");
        // Create a list of rmax values to test
        List<double> rmaxValues = new List<double>();
        for (double rmax = 5.0; rmax <= 50.0; rmax += 5.0)
        {
            rmaxValues.Add(rmax);
        }
        // Calculate ground state energy for each rmax
        List<double> energies = new List<double>();
        foreach (double rmax in rmaxValues)
        {
            double energy = CalculateGroundStateEnergy(rmax, fixedDr);
            energies.Add(energy);
            Console.WriteLine($"rmax = {rmax}, ε₀ = {energy}");
        }
        // Write results to a file for plotting
        using (StreamWriter writer = new StreamWriter("rmax_convergence.dat"))
        {
            writer.WriteLine("# rmax ε₀");
            for (int i = 0; i < rmaxValues.Count; i++)
            {
                writer.WriteLine($"{rmaxValues[i]} {energies[i]}");
            }
        }
        Console.WriteLine("Results saved to rmax_convergence.dat");
    }
 static double CalculateGroundStateEnergy(double rmax, double dr){
		int npoints =(int)(rmax/dr)-1;
		vector r = new vector(npoints);
		for(int i=0;i<npoints;i++)r[i]=dr*(i+1);
		matrix H = new matrix(npoints,npoints);
		for(int i=0;i<npoints-1;i++){
   			H[i,i]  =-2*(-0.5/dr/dr);
   			H[i,i+1]= 1*(-0.5/dr/dr);
   			H[i+1,i]= 1*(-0.5/dr/dr);
  			}
		H[npoints-1,npoints-1]=-2*(-0.5/dr/dr);
		for(int i=0;i<npoints;i++)H[i,i]+=-1/r[i];
		var HdiagwV = jacobi.cyclic(H);
		var Evals = HdiagwV.Item2;
		var Evecs = HdiagwV.Item3;
		double eps_0=Evals[0];
		return eps_0;
		}
	}

