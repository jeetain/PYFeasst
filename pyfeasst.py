import os, sys

def ntot_window_scaling (n_f, dw, w_max, n_ov):
	"""
	Compute the upper bounds for windows of a flat histogram simulation using N_{tot} as the order parameter.
	
	Parameters
	----------
	n_f : int
		Max total number of particles (n_tot)
	dw : int
		Final window width
	w_max : int
		Number of windows to use
	n_ov : int
		Number of overlapping points to use
    
	Returns
	-------
	ndarray
		Array of tuples of (lower, upper) bound for each window
   
	"""

	dw -= n_ov # account for overlap
	assert (n_ov < w_max), "n_ov too large"

	alpha = np.log(float(n_f)/(float(n_f) - float(dw))) / np.log(w_max/(w_max-1.0))
	coeff = float(n_f)/(float(w_max)**alpha)

    
	x = np.linspace(1, w_max, w_max)
	ub = np.round(coeff*x**alpha).astype(int)
	lb = [0]
	for i in xrange(1, int(w_max)):
		lb.append(ub[i-1]-n_ov+1)

	return zip(lb, ub)


class Base (object):
	def __init__ (self):
		"""
		Initialize the class by importing FEASST.

		"""

                global feasst
                self.get_feasst()
                import feasst

	def get_feasst (self):
		"""
		Find the FEASST install directory and add it to the python path.
		This requiers FEASST_INSTALL_DIR_ to be defined in the shell environment.

		"""

	        try:
        	        feassthead = os.getenv("FEASST_INSTALL_DIR_")
        	except:
                	print 'Unable to locate FEASST install directory'
                	sys.exit()

        	def find_all(name, path):
                	result = []
                	for root, dirs, files in os.walk(path):
                        	if name in files:
                                	result.append(os.path.join(root, name))
                	return result

        	feasstname = find_all("_feasst.so", feassthead)
        	if (len(feasstname) > 1):
                	print ('Found multiple FEASST sources which may cause an error: ', feasstname)
                	sys.exit()
        	else:
                	feasstdir = feasstname[0].split('/_feasst.so')[0]
                	sys.path.append(feasstdir)

	def run (self):
		"""
		Run a simulation.

		"""

		print "Beginning simulation"

	def restart (self):
		"""
		Restart a simulation.

		"""

		print "Restarting simulation"
