using static System.Console;
using static System.Math;
public class Program{
static int Main(){
	var rnd=new System.Random();
	var u=new vec(rnd.NextDouble(),rnd.NextDouble(),rnd.NextDouble());
	var v=new vec(rnd.NextDouble(),rnd.NextDouble(),rnd.NextDouble());
	u.print("u=");
	v.print("u=");
	WriteLine($"u={u}");
	WriteLine($"v={v}");
	WriteLine();
	vec t;

	t=new vec(-u.x,-u.y,-u.z);
	(-u).print("-u =");
	t.print   ("t  =");
	if(vec.approx(t,-u))WriteLine("test 'unary -' passed\n");

	t=new vec(u.x-v.x,u.y-v.y,u.z-v.z);
	(u-v).print("u-v =");
	t.print    ("t   =");
	if(vec.approx(t,u-v))WriteLine("test 'operator-' passed\n");

	t=new vec(u.x+v.x,u.y+v.y,u.z+v.z);
	(u+v).print("u+v =");
	t.print    ("t   =");
	if(vec.approx(t,u+v))WriteLine("test 'operator+' passed\n");

	double c=rnd.NextDouble();
	t=new vec(u.x*c,u.y*c,u.z*c);
	var tmp=u*c; // bug in mcs
	tmp.print("u*c =");
	t.print  ("t   =");
	if(vec.approx(t,u*c))WriteLine("test 'operator*' passed\n");

	double d = u.x*v.x+u.y*v.y+u.z*v.z;
	double dtest =vec.dot(u,v);
	WriteLine($"u dot v= {dtest}");
	WriteLine($"d  = {d}");
	if( vec.approx(d, dtest) )WriteLine("test 'operator%' passed\n");

	return 0;
	}//main
}//program
