#f(x) = sin(3.14159*(x+0.5*(2*x**2-x**3)))
#g(x) = sin(3.14159*(x+0.5*(2*x**2-x**3)))
f(x) = sin(3.14159*(x+0.5*(x**2-x**3)))
g(x) = sin(3.14159*(x+0.5*(x**2-1.1*x**3)))
#h(x) = sin(3.14159*(x+0.88*(4*x**2-x**4)))
#i(x) = sin(3.14159*(x+0.89*(4*x**2-x**4)))
h(x) = sin(3.14159*(x+0.03*(8*x-x**4)))
i(x) = sin(3.14159*(x+0.02*(8*x-x**4)))
h(x) = sin(3.14159*(x+0.01*(8*x-x**4)))
i(x) = sin(3.14159*(x+0.00*(8*x-x**4)))
plot [0:2] f(x), g(x), h(x), i(x)
