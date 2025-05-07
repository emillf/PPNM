using static System.Math;
public static class sfunsln{
public static double lnfgamma(double x){ ///single precision gamma function
        if(x <= 0)return double.NaN
        if(x<9)return lnfgamma(x+1)/x; // Recurrence relation
        double lnfgamma=x*Log(x+1/(12*x-1/x/10))-x+Log(2*PI/x)/2;
        return lnfgamma;
        }
}
