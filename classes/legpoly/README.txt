README.txt - Class "legpoly"

A (probably not so efficient) class to model Legendre polynomials, to be used 
with QUINOPT.

The actual class is put as a subfolder because the functions legMulx and legProd
are not methods for the class, and must accept arguments of other classes. 

However, these functions are specific to the "legpoly" class, and are only used by
methods of the "legpoly" class.