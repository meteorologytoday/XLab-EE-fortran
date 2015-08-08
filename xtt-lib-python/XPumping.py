import numpy as np;

class Pumping:
	"""
	Here "w" actually means "rho w", but rho is not important during the
	calculation process, so it is free of users to treat "w" in their own
	ways.
	
	w(0 < r < r_arr[0]) = 0
	w(r_arr[0] < r < r_arr[1]) = coe[0][0] * (r - r_arr[0]) * (r - r_arr[1])
	w(r_arr[1] < r < r_arr[2]) = coe[1][0] * (r - r_arr[1]) * (r - r_arr[2])
	w(r > r_arr[2]) = 0
	
	coe[0][0] = - 4 * w0 / (r_arr[1] - r_arr[0])**2
	
	# coe[1][0] is contrained by coe[0][0]
	
	psi(0 < r < r_arr[0]) = 0
	psi(r_arr[0] < r < r_arr[1]) = coe[0][1] + coe[0][0] * int_part(r, r_arr[0], r_arr[1]);
	psi(r_arr[1] < r < r_arr[2]) = coe[1][1] + coe[1][0] * int_part(r, r_arr[1], r_arr[2]);
	psi(r > r_arr[2]) = 0
	
	# constrain
	(a) psi(r_arr[0]) continuous
	(b) psi(r_arr[1]) continuous
	(c) psi(r_arr[2]) continuous
	
	coe[0][1] = - coe[0][0] * int_part(r_arr[0], r_arr[0], r_arr[1]);
	coe[1][0] * int_part(r_arr[2], r_arr[1], r_arr[2]) + coe[1][1] = 0;
	coe[1][0] * int_part(r_arr[1], r_arr[1], r_arr[2]) + coe[1][1] = coe[1][0] * int_part(r_arr[1], r_arr[0], r_arr[1]) + coe[0][1];
	"""
	def __init__(self, rho_w0, r_arr):
		if(len(r_arr) != 3):
			raise Exception("The length of r array must be exactly 3, the input length is %d" % (len(r_arr),));

		self.rho_w0 = rho_w0;
		self.r_arr = r_arr.copy();
		self.init_const();

	def int_part(self, at_r, r_min, r_max):
		return ((at_r**4.0) / 4.0 - (r_min + r_max) / 3.0 * (at_r**3.0) + r_min * r_max * (at_r ** 2.0) / 2.0 );
	
	
	
	def int_part_up(self, at_r):
		return self.int_part(at_r, self.r_arr[0], self.r_arr[1]);

	def getTotalFlux(self):
		return self.coe[0][0] * (self.int_part_up(self.r_arr[1]) - self.int_part_up(self.r_arr[0]));
		
	def getFluxGeometry(self):
		""" Geometry factor, which means flux_geometry * w0 = total upward flux (positive) """
		return self.getTotalFlux() / self.rho_w0;
	

	
		
	def init_const(self):
		self.coe = np.array([[0.0,0.0],[0.0,0.0]]);
		
		self.coe[0][0] = - 4.0 * self.rho_w0 / (self.r_arr[1] - self.r_arr[0])**2.0;
		self.coe[0][1] = - self.coe[0][0] * self.int_part(self.r_arr[0], self.r_arr[0], self.r_arr[1]);
		
		# ax = b
		a = np.array([ 
			[self.int_part(self.r_arr[2], self.r_arr[1], self.r_arr[2]), 1],
			[self.int_part(self.r_arr[1], self.r_arr[1], self.r_arr[2]), 1]
		]);
		
		b = np.array([
			0,
			self.coe[0][0] * self.int_part(self.r_arr[1], self.r_arr[0], self.r_arr[1]) + self.coe[0][1]
		]);
		
		x = np.linalg.solve(a, b);
		self.coe[1][0], self.coe[1][1] = x[0], x[1];
		
	
	def getRPsi(self, r):
		rPsi = 0;
		if(r <= self.r_arr[0]):
			rPsi = 0;
		elif(r <= self.r_arr[1]):
			rPsi = self.coe[0][0] * self.int_part(r, self.r_arr[0], self.r_arr[1]) + self.coe[0][1];
		elif(r <= self.r_arr[2]):
			rPsi = self.coe[1][0] * self.int_part(r, self.r_arr[1], self.r_arr[2]) + self.coe[1][1];
		else:
			rPsi = 0;
			
		return rPsi;
	
	def getRhoW(self, r):
		rho_w = 0;
		if(r <= self.r_arr[0]):
			rho_w = 0;
		elif(r <= self.r_arr[1]):
			rho_w = self.coe[0][0] * (r - self.r_arr[0]) * (r - self.r_arr[1]);
		elif(r <= self.r_arr[2]):
			rho_w = self.coe[1][0] * (r - self.r_arr[1]) * (r - self.r_arr[2]);
		else:
			rho_w = 0;
			
		return rho_w;