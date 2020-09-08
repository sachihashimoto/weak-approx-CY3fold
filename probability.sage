
def size_Gr(k, n, q):
	'''computes the number of F_q points on the
	Grassmannian G(k, n) of *projective* k-planes
	in projective n-space over F_q'''
	return q_analog_factorial(n+1, q)/(q_analog_factorial(k+1, q)*q_analog_factorial(n-k, q))

def q_analog(n, q):
	return (q^n - 1)/(q-1)

def q_analog_factorial(n, q):
	F = 1
	for k in range(1, n+1):
		F *= q_analog(k, q)

	return F


def bad_set_p(p):
	'''computes the sum of the number of rank 1 5x5 matrices
	and the number of rank 2 nonsplit 5x5 matrices over F_p 
	up to common scaling, WHEN p IS ODD'''

	num_rk1 = size_Gr(3, 4, p) 	#in rank 1 we just have choice of the kernel

	num_rk2 = size_Gr(2, 4, p) * (p-1)/2 * p 	#in rank 2, we have the kernel, discriminant, and one more free param

	return num_rk1 + num_rk2

def checkSmoothPoints(Uc,Ucdx,x):
	'''INPUT: a quadric Uc over char 2 in five variables x, and list of derivatives Ucdx, and the list of variables x
	OUTPUT: True if smooth point over F2, FALSE otherwise
	'''
	for pt in product(range(2), repeat = 5):
		substitutions = dict(zip(x,pt))
		Uceval = Uc.subs(substitutions)
		if Uceval == 0:
			for Ucdxi in Ucdx:
				if Ucdxi.subs(substitutions) != 0:
					return True
	return False

def bad_set_2():
	'''	computes the number of quadrics in 5 variables 
	in characteristic 2 that do not have a smooth F_2 point
	this function takes several minutes to return a value
	'''
	R.<x0,x1,x2,x3,x4,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14>=GF(2)[]
	x = [x0,x1,x2,x3,x4]
	c = [c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14]
	BadSet = 0
	U = sum(x*coeff for (x,coeff) in zip([x0^2 , x1^2 , x2^2 , x3^2 ,x4^2 , x0*x1, x0*x2,x0*x3,x0*x4,x1*x2,x1*x3,x1*x4,x2*x3,x2*x4,x3*x4],c))
	for coeffs in product(range(2),repeat = 15):								
		if not any(coeffs):
			continue
		substitutions = dict(zip(c,coeffs))
		Uc = U.subs(substitutions)
		Ucdx = [ Uc.derivative(xi) for xi in x]
		if not any(Ucdx):
			BadSet = BadSet+1
			continue
		smPt = checkSmoothPoints(Uc,Ucdx,x)
		if not smPt:
			BadSet = BadSet+1
	return BadSet

def bad_factor_at_p(p):
	'''Computes the bad factor at p for p > 2
	as described in Prop 3.4 of Hashimoto, Honigs, Lamarche, Vogt'''
	if p == 2:
		return (bad_set_2()*size_Gr(3, 13, 2))/size_Gr(4, 14, 2)

	else:
		return (bad_set_p(p)*size_Gr(3, 13, p))/size_Gr(4, 14, p)

def factor_at_p(p):

	return 1 - bad_factor_at_p(p)


#we can compute the bound in Proposition 3.4
prod = 1
for p in prime_range(100):
	prod = prod * factor_at_p(p)
print(prod)


