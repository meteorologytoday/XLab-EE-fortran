class WindProfile:
	"""WindProfile needs Coriolis parameter and its boundary radius, which should be
	one element less then previous"""
	def __init__(self, f0, f_arr, radius_arr):
		self.f0 = f0;
		self.f_arr = f_arr[:];
		self.radius_arr = radius_arr[:];
		self.makeProfile();
	
	def makeProfile(self):
		self.konst = [0.0 for _dummy in self.f_arr];
		for i in range(1, len(self.konst)):
			self.konst[i] = self.konst[i-1] + (self.radius_arr[i-1]**4.0) / 4.0 * (self.f_arr[i-1]**2.0 - self.f_arr[i]**2.0);
			#print("konst[%d] = %f" % (i,self.konst[i]));
	
	def getWind(self, r):
		region = len(self.f_arr) - 1;
		for i in range(0, len(self.radius_arr)):
			if(r < self.radius_arr[i]):
				region = i;
				break;
		
		return (r**2.0 / 4.0 * self.f_arr[region]**2.0 + self.konst[region]/r**2.0)**(1.0/2.0) - 1.0/2.0 * self.f0 * r if(r != 0.0) else 0.0;
