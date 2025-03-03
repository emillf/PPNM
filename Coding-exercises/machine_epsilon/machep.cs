using System;
public class Machep{
	public double epvaldouble{ get; private set; }
	public float epvalfloat{ get; private set; }
	public double epvalaprxdouble{ get; private set; }
	public float epvalaprxfloat{ get; private set; }
	public double tiny{ get; private set; }
	public double a{ get; private set; }
	public double b{ get; private set; }
	public Machep()
	{
		double x=1;
		while(1+x!=1)
		{
		x/=2;
		}
		x*=2;
		float y=1F;
		while((float)(1F+y) != 1F)
		{
		y/=2F;
		}
		y*=2F;
		epvaldouble = x;
        	epvalfloat = y;
        	epvalaprxdouble = Math.Pow(2,-52);
        	epvalaprxfloat = (float)Math.Pow(2,-23);
		tiny = Math.Pow(2,-52)/2;
		a=1+tiny+tiny;
		b=tiny+tiny+1;
	}
}
