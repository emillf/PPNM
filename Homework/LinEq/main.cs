using System;
using static System.Console;
public class Program{
	static void Main(){
        	var rnd = new Random(3);
        	double x = rnd.NextDouble();
        	double y = rnd.NextDouble();
        	double z = rnd.NextDouble();
		var A = matrix.random(7,3,rnd);
		A.print("Random matrix A");
		WriteLine("/n A decomposed into Q and R gives:/n");
		QR = new QR;
		var QRfact = QR.decomp(A);
		var Q = QRfact.Item1;
		var R = QRfact.Item2;
		WriteLine("Q = /n");
		Q.print;
		WriteLine("R = /n");
		R.print;
		bool AeqQR=A.approx(QRfact);
		WriteLine($"QR = A? {AeqQR}/n");
		}
	}
