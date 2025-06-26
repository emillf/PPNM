
NOTE: The program does take a while to run, so please wait at least 10-20 seconds after make clean && make

ALSO READ Out.txt

Part A)
To do part A required creating a method to create an analytic Jacobian on the form suggested in the exercise.
This was done simply following the form given, and taking the final column of J to contain the derivatives of λ.

Secondly, to easily adapt my earlier linesearch method I needed to create a method that output the function we wished to minimize,
namely f(x)={Av-λv,v^Tv-1} where x=(v,λ). This was done thorugh Create_func_from_A taking only the matrix A and outputting f.
 
Finally, I adapted my newton linesearch to instead take my startguesses for v and λ, as well as the symmetric matrix A.
This meant that i could constuct f as needed. I of course also needed to collect my startguesses into one vector x to match the Jacobian
and f.

To test, I generated a random 6x6 symmetric matrix A, and saw that for the eigenvectors and eigenvalues found Av=λv. 

Part B)
For this part I simply needed to add a stopwatch at the beginning and at the end of my Newton eigenvaluefinder. 
I then timed the procedure for dimensions up to nxn=400 in large increments, as to not spend too much computation time. After n=400
the computation time was too long.

To get a bit more consistency, and to closer match the later H, I chose a Toeplits tridiagonal matrix to test on.
I plotted t as a function of n on a log scale to see that it did indeed seem to go like O(n^2).

Part C)
This part first required constructing a numerical Hamiltonian with which to find eigenvalues and eigenvectors of.

I constructed this H by approximating the second derivative in the radial Schrödinger equation for hydrogen with reduced radial 
wavefunctions using the central difference approximation. This yielded terms of H containing u(r) and u(r±2δr). The terms containing 
u(r±2δr) would then be the off diagonal terms of H while the diagonal terms contained u(r). I also used an option for creating 
a different H depending on l. 

Secondly, due to both the hint given, and issues with divergence in the numerical method, I tried implementing an option for
analytical linesearch. This was done by Taylor expanding f(x-αΔx) to first order. f(x-αΔx)≈f(x)+αJΔx. Then writing 
||f(x)+αJΔx||^2=(f(x)+αJΔx)^T*(f(x)+αJΔx) and solving for df/dα=0, I get the formula α=-f^T*J*Δx/||JΔx||^2. Thus i find the alpha that
minimizes ||f(x)+αJΔx||^2.

Finally I do the linesearch to find eigenvectors and eigenvalues of H_(l=0) and H(l=1). I chose initial energies to be 1/(2n^2+0.1) and
u_i - ureal_i = 0.1. I could not find a way to implement H such as to get better eigenvalues of the correct behaviour of u.

 
Self-evaluation:
Part A and Part B were completed without any issues, so that would be 9 points. Part C) was not sufficiently completed in my opinion.
I did investigate and implement an analytical linesearch, create a somewhat correct Hamiltonian, and plot several lowest eigenfunctions
and find their energies. But the energies were only somewhat correct, and the eigenfunctions that I reached did not display the desired 
behaviour close to r=0. I think this would net me 0.5 points, as I consider it a task half completed.

Final score 9.5/10.  
