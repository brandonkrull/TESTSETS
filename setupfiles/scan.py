"""
Called by ./slurm.dscf.realscan
"""
import shutil as sh
import subprocess as sp
import re
import numpy as np

# print "inside python"
# f = open("control", "r")
# text = f.read()
# text = re.sub("\$scfiterlimit.*", "$scfiterlimit 1", text)
# f.close()
# f = open("control", "w")
# f.write("$acgga ")
# f.write(text)
# f.close()

for fun in ['acgga0', 'b-acgga0']:
    for exx in np.arange(0.12, 0.23, 0.01):
        try:
            sh.copy("ini/mos", ".")
        except IOError:
            print "Unrestricted calculation. No mos"
            sh.copy("ini/alpha", ".")
            sh.copy("ini/beta", ".")

        f = open("control", "r")
        text = f.read()
        text = re.sub("functional.*", "functional " + fun, text)
        f.close()
        f = open("control", "w")
        f.write(text)
        f.close()

        f = open("control", "r")
        text = f.read()
        text = re.sub("\$acgga.*", "$acgga " + str(exx), text)
        f.close()
        f = open("control", "w")
        f.write(text)
        f.close()

        suffix = "{:.3f}".format(exx).split(".")[1]
        out  = sp.check_output(["dscf"])
        fout = open("dscf." + fun + "." + suffix, "w")
        fout.write(out)
        fout.close()
