using static System.Math;
public static class sfuns{
public static double fgammaln(double x){ ///single precision gamma function
        if(x <= 0)return double.NaN;
        if(x<9)return fgammaln(x+1)/x; // Recurrence relation
        double lnfgamma=x*Log(x+1/(12*x-1/x/10))-x+Log(2*PI/x)/2;
        return lnfgamma;
        }
}
