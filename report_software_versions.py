# -*- coding: utf-8 -*-
import sys, importlib
print("Python", sys.version.split()[0])
for m in ["numpy", "scipy", "pandas", "matplotlib", "h5py", "mofapy2", "openpyxl"]:
    try:
        mod = importlib.import_module(m)
        print(m + "==" + str(getattr(mod, "__version__", "unknown")))
    except Exception:
        print(m + " NOT INSTALLED")
