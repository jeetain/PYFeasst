import os, sys
sys.path.append("../../py/")
import cylinder_fslj as cfslj

if __name__ == "__main__":
	x = cfslj.cylinder_fslj()
	x.run(steps=int(1e9), orderMin=215, orderMax=280, temp=1.35, mu1=0.0, dmu2=1.65, order=3)
