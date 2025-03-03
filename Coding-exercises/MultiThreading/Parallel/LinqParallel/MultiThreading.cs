using System;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using static System.Console;

public class Program
{
    public static void Main(string[] args)
    {
        int nterms = (int)1e8; // default value

        // Parse command-line arguments
        foreach (var arg in args)
        {
            var words = arg.Split(':');
            if (words[0] == "-terms") nterms = (int)float.Parse(words[1]);
        }

        // Thread-local storage for partial sums
        var sum = new ThreadLocal<double>(() => 0.0, trackAllValues: true);

        // Parallel computation of the harmonic series
        Parallel.For(1, nterms + 1, i =>
        {
            sum.Value += 1.0 / i;
        });

        // Sum all the partial sums from each thread
        double totalsum = sum.Values.Sum();

        // Output the result
        WriteLine($"Harmonic sum up to {nterms} is {totalsum} using automatic local linear quadratic threading");
    }
}
