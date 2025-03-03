using System;
using System.Threading;
using static System.Console;

public class data { public int a,b; public double sum;}
public class Program{
	public static void harm(object obj){
        	var arg = (data)obj;
        	arg.sum=0;
        	for(int i=arg.a;i<arg.b;i++) arg.sum+=1.0/i;
        	}
	public static void Main(string[] args){
		int nthreads = 1, nterms = (int)1e8; /* default values */
                        foreach(var arg in args) {
                        var words = arg.Split(':');
                        if(words[0]=="-threads") nthreads=int.Parse(words[1]);
                        if(words[0]=="-terms"  ) nterms  =(int)float.Parse(words[1]);
                        }
		data[] Params = new data[nthreads];
		for(int i = 0 ; i<nthreads ; i++){
			Params[i] = new data();
			Params[i].a = (nterms/nthreads)*(i)+1;
			Params[i].b = ((nterms)/nthreads)*(i+1);
			}
		Params[nthreads-1].b=nterms+1;
		var threads = new Thread[nthreads];
		for(int i = 0 ; i<nthreads ; i++){
			threads[i] = new Thread(harm);
			threads[i].Start(Params[i]);
			}
		foreach(var thread in threads) thread.Join();
		double total=0; foreach(var p in Params) total+=p.sum;
		WriteLine($"Harmonic sum up to {nterms} is {total} using {nthreads} threads");
		}
	}
