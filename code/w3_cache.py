#!/usr/bin/env python3
"""Run tf_path_check.py (unchanged, 21 eps points plus its refinement) and cache Z_(N,l)(eps) to w3_path.npz.
Takes about 5 to 10 minutes.  nice -n 19 env OMP_NUM_THREADS=1 python3 w3_cache.py"""
import os, sys, runpy
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = HERE
sys.argv = ["tf_path_check.py", "21"]
sys.path.insert(0, CODE)
g = runpy.run_path(os.path.join(CODE, "tf_path_check.py"), run_name="w3")
rows = [(e, N, l, v[0]) for e, zd in g["Zs"].items() for (N, l), v in zd.items()]
np.savez(os.path.join(HERE, "w3_path.npz"), rows=np.array(rows), a=g["a"], A=g["A"])
print("saved", len(rows))
