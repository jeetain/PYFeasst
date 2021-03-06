import os, sys, math
sys.path.append("${HOME}/PYFeasst/")
import pyfeasst as pf

class cylinder_fslj (pf.Base):
	"""
	Binary linear force-shifted Lennard-Jones mixture in cylindrical confinment (in X-direction).

	"""

	def __init__ (self):
		"""
		Initialize the class.

		"""

		global feasst
		self.get_feasst()
		import feasst

	def run (self, orderMax=100, orderMin=0, temp=1.50, mu1=0.0, dmu2=0.0, order=3, steps=10000):
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

		"""

		s = feasst.Space(3,0)
		s.lset(15)
		p = feasst.PairLJMulti(s, s.minl()/4.0)

		molA = os.path.dirname(os.path.realpath(__file__))+"/data.lj"
		molB = os.path.dirname(os.path.realpath(__file__))+"/data.ljb"

		s.addMolInit(molA)
		p.initLMPData(molA)
		s.addMolInit(molB)
		p.initLMPData(molB)

		p.rCutijset(0,0,(p.sig(0)+p.sig(0))/2*3.0)
		p.rCutijset(0,1,(p.sig(0)+p.sig(1))/2*3.0)
		p.rCutijset(1,1,(p.sig(1)+p.sig(1))/2*3.0)

		for i in range(2):
                        for j in range(2):
                                p.linearShiftijset(i,j,True)

		s.updateCells(p.rCutMaxAll()) # intermolecular potential needs max rcut

		barrier = feasst.SpeciesBarriers (s.nParticleTypes())
		radius = 4.0 # this goes to CENTER OF MASS, so pore radius actually = 4.5 with sigma/2 exclusion
		eps_w = 3.0
		width = 0.5
		pt = [0.0, 0.0, 0.0]
		barrier.addSqwCylinder (0, pt, radius, eps_w, width, 0)
		barrier.addSqwCylinder (1, pt, radius-(p.sig(1)-p.sig(0))/2, eps_w/2.0, width, 0)
		pairBarr = feasst.PairBarriers (s, barrier)

		ipair = feasst.PairHybrid (s, max(p.rCutMaxAll(), width))
		ipair.addPair(p)
		ipair.addPair(pairBarr)
		ipair.initEnergy()

		assert (barrier.safe(ipair)) # check barrier doesn't violate any boundary conditions or allow periodic interactions

		orderParam = "nmol"
		beta = 1.0/temp

		cc = feasst.CriteriaWLTMMC (beta, math.exp(beta*mu1), orderParam, orderMin-0.5*1.0, orderMax+0.5*1.0, orderMax-orderMin+1)
		cc.addActivity(math.exp(beta*(mu1+dmu2)))
		mc = feasst.WLTMMC (s, ipair, cc)
		em = feasst.AnalyzeExtensiveMoments(s, ipair)
		em.setOrderParam (order, cc)

		mc.weight = 0.6
		feasst.deleteTrial(mc, molA)
		mc.weight = 0.6
		feasst.deleteTrial(mc, molB)
		mc.weight = 0.6
		feasst.addTrial(mc, molA)
		mc.weight = 0.6
		feasst.addTrial(mc, molB)
		mc.weight = 0.2
		feasst.xswapTrial(mc)
		mc.weight = 0.4
		feasst.transformTrial(mc, "translate", 0.4)

		mc.nMolSeek(orderMin, molA, 1000000) # mininum bound with pure A (smaller of the 2 species)

		nfreq = 100000
		mc.initLog("log", nfreq)
		mc.initColMat("colMat", nfreq)
		mc.setNFreqCheckE(nfreq, 1e-8)
		mc.setNFreqTune(nfreq)
		mc.initMovie("movie", nfreq*1000)
		mc.initRestart("tmp/rst", nfreq*1000)

		em.initFreq(5)
		em.initPrintFreq(nfreq*1000)
		em.initFileName("extMom")
		mc.initAnalyze(em)

		cc.collectInit(4.0e-6) # start collecting TMMC matric at wlFlat == 18
		cc.tmmcInit(5.0e-7)  # switch to TMMC after 22
		mc.wlFlatProduction(22) # after this # iterations, switch to "production" phase and measure
		mc.runNumTrials(steps)

		em.write()

if __name__ == "__main__":
	x = cylinder_fslj()
	x.run(steps=0)
