import os, sys
sys.path.append("../../py/")
import binary_fslj as bfslj

if __name__ == "__main__":
	x = bfslj.binary_fslj()
	x.run(steps=int(1e7), orderMin=319, orderMax=393, temp=1.50, mu1=0.0, dmu2=1.65, order=3)
