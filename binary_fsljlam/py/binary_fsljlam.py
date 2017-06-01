import os, sys, math
sys.path.append("${HOME}/PYFeasst/")
import pyfeasst as pf

class binary_fsljlam (pf.Base):
	"""
	Binary linear force-shifted Lennard-Jones Lambda mixture using LB mixing rules.

	"""

	def __init__ (self):
		"""
		Initialize the class.

		"""

		global feasst
		self.get_feasst()
		import feasst

	def parse (**kwargs):
		self.orderMax = kwargs.get("orderMax", 100)
		self.orderMin = kwargs.get("orderMax", 0)
		self.temp = kwargs.get("temp", 1.50)
		self.mu1 = kwargs.get("mu1", 0.0)
		self.dmu2 = kwargs.get("dmu2", 0.0)
		self.order = kwargs.get("order", 3)
		self.steps = kwargs.get("steps", 10000)
		self.lam11 = kwargs.get("lam11", 0.0)
		self.lam22 = kwargs.get("lam22", 0.0)
		self.molA = kwargs.get("molA", os.path.dirname(os.path.realpath(__file__))+"/data.lj")
		self.molB = kwargs.get("molB", os.path.dirname(os.path.realpath(__file__))+"/data.ljb")

	def run (self, **kwargs):
		"""
		Run a new simulation from scratch.

		Parameters
		----------
		orderMax : int
			Max number of particles to allow in this window (default=100)
		orderMin : int
			Min number of particles to allow in this window (default=0)
		temp : double
			Temperature (default = 1.5)
		mu1 : double
			Chemical potential of species 1
		dmu2 : double
			mu2 - mu1
		order : int
			Max order to record extensive moments up to
		steps : int
			Number of MC steps to perform
		lam11 : float
			Value of lambda_11
		lam22 : float
			Value of lambda_22
		molA : string
			Filename for molA (1)
		molB : string
			Filename for molB (2)

		"""

		self.parse(kwargs)

		s = feasst.Space(3,0)
		s.lset(9.0)
		p = feasst.PairLJMulti(s, s.minl()/2.0)

		s.addMolInit(self.molA)
		p.initLMPData(self.molA)
		s.addMolInit(self.molB)
		p.initLMPData(self.molB)

		# Force-shift potential for each pair
		p.rCutijset(0,0,(p.sig(0)+p.sig(0))/2*3.0)
		p.rCutijset(0,1,(p.sig(0)+p.sig(1))/2*3.0)
		p.rCutijset(1,1,(p.sig(1)+p.sig(1))/2*3.0)

		p.setLambdaij(0,0, self.lam11)
		p.setLambdaij(0,1, 0.5*(self.lamm11+self.lam22)) # LB mixing
		p.setLambdaij(1,1, self.lam22)
		p.lrcFlag = 0;

		for i in range(2):
			for j in range(2):
				p.linearShiftijset(i,j,True)

		s.updateCells(p.rCutMaxAll()) # intermolecular potential needs max rcut

		ipair = feasst.PairHybrid (s, p.rCutMaxAll())
		ipair.addPair(p)
		ipair.initEnergy()

		orderParam = "nmol"
		beta = 1.0/self.temp

		cc = feasst.CriteriaWLTMMC (beta, math.exp(beta*self.mu1), orderParam, self.orderMin-0.5*1.0, self.orderMax+0.5*1.0, self.orderMax-self.orderMin+1)
		cc.addActivity(math.exp(beta*(self.mu1+self.dmu2)))
		mc = feasst.WLTMMC (s, ipair, cc)
		em = feasst.AnalyzeExtensiveMoments(s, ipair)
		em.setOrderParam (self.order, cc)

		mc.weight = 0.6
		feasst.deleteTrial(mc, self.molA)
		mc.weight = 0.6
		feasst.deleteTrial(mc, self.molB)
		mc.weight = 0.6
		feasst.addTrial(mc, self.molA)
		mc.weight = 0.6
		feasst.addTrial(mc, self.molB)
		mc.weight = 0.2
		feasst.xswapTrial(mc)
		mc.weight = 0.4
		feasst.transformTrial(mc, "translate", 0.4)

		mc.nMolSeek(self.orderMin, self.molA, 1000000) # mininum bound with pure A (smaller of the 2 species)

		nfreq = 100000
		mc.initLog("log", nfreq)
		mc.initColMat("colMat", nfreq*1000)
		mc.setNFreqCheckE(nfreq*1000, 1e-8)
		mc.setNFreqTune(nfreq*1000)
		mc.initMovie("movie", nfreq*1000)
		mc.initRestart("tmp/rst", nfreq*1000)

		em.initFreq(5) # from 1 -> 5 makes this as efficient as measuring 2nd vs 3rd order moments
		em.initPrintFreq(nfreq*1000)
		em.initFileName("extMom")
		mc.initAnalyze(em)

		cc.collectInit(4.0e-6) # start collecting TMMC matric at wlFlat == 18
		cc.tmmcInit(5.0e-7)  # switch to TMMC after 22
		mc.wlFlatProduction(22) # after this number of iterations, switch to "production" phase and measure
		mc.runNumTrials(self.steps)

		em.write() # write at the end to ensure up-to-date

if __name__ == "__main__":
	x = binary_fsljlam()
	x.run(steps=0)
