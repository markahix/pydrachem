import imp
libraries=['glob','matplotlib','numpy','scipy','prody']
for lib in libraries:
    try:
        imp.find_module(lib)
    except:
        print(lib,"dependency missing.")
   
try:
    for lib in libraries:
        imp.find_module(lib)
    from .correlation import *
    from .rmsf import *
    from .rmsd import *
    from .normalmodes import *
    from .hbonds import *
    from .eda import *
except:
    print("Dependencies missing.  'Subplots' submodule not loaded.")