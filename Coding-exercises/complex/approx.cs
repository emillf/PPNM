namespace MathLibrary
{
	public static class MathUtils
	{
		public static bool Approx
		(double a, double b, double acc=1e-9, double eps=1e-9)
		{
			if(System.Math.Abs(b-a) <= acc) return true;
			if(System.Math.Abs(b-a) <= System.Math.Max(System.Math.Abs(a),System.Math.Abs(b))*eps) return true;
			return false;
		}
	}
}
