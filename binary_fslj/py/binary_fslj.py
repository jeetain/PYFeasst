import os, sys, math
sys.path.append("${HOME}/PYFeasst/")
import pyfeasst as pf

class binary_fslj (pf.Base):
	"""
	Binary linear force-shifted Lennard-Jones mixture.

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
			Min number of particles to allow in this window (default=100)
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
		s.lset(10)
		p = feasst.PairLJMulti(s, s.minl()/2.0)

		molA = os.path.dirname(os.path.realpath(__file__))+"/data.lj"
		molB = os.path.dirname(os.path.realpath(__file__))+"/data.ljb"

		s.addMolInit(molA)
		p.initLMPData(molA)
		s.addMolInit(molB)
		p.initLMPData(molB)

		p.rCutijset(0,0,(p.sig(0)+p.sig(0))/2*3.0)
		p.rCutijset(0,1,(p.sig(0)+p.sig(1))/2*3.0)
		p.rCutijset(1,1,(p.sig(1)+p.sig(1))/2*3.0)
		s.updateCells(p.rCutMaxAll()) # intermolecular potential needs max rcut

		p.linearShift(True) # force-shifted at rcut so dU/dr = 0
		p.cutShift(True) # also shift U = 0 at rcut

		ipair = feasst.PairHybrid (s, p.rCutMaxAll())
		ipair.addPair(p)
		ipair.initEnergy()

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

		em.initFreq(5) # from 1 -> 5 makes this as efficient as measureing 2nd vs 3rd order moments
		em.initPrintFreq(nfreq*1000)
		em.initFileName("extMom")
		mc.initAnalyze(em)

		cc.collectInit(4.0e-6) # start collecting TMMC matric at wlFlat == 18
		cc.tmmcInit(5.0e-7)  # switch to TMMC after 22
		mc.wlFlatProduction(22) # after this # iterations, switch to "production" phase and measure
		mc.runNumTrials(steps)

		em.write() # write at the end to ensure up-to-date

if __name__ == "__main__":
	x = binary_fslj()
	x.run(steps=0)
