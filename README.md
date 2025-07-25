# POMMES v1.1
## Peak-Ordering Measurement for Multi-Elemental Systems

POMMES is a tool written in Python 3 for quantifying order in multi-elemental crystal systems such as perovskite crystals. We aim to deliver a GUI capable of estimating order in multiple crystal sites given ratios between peak intensities in data from X-ray diffraction (XRD).

`main.py` plots peak ratios versus order from two .cif files, one for the ordered and one for the disordered state. If given experimentally determined peak ratio, it can interpolate them to estimate the order in the sample.
`GUI.py` does the same task as above but with a helpful GUI.
`xrdmlplotter.py` can plot XRD data from .xrdml or .csv files and extract peak data, and users can use it to calculate peak ratios required in `main.py`.
`readCIF.py`, if run on two .cif files as above, will output a snippet of code that can be pasted into `preload.py` representing the necessary data from both .cifs. By setting `preload_matrices=True` in `main.py`, data will be loaded from `preload.py` instead of rereading the .cifs every time the code is run, saving time and computational power.

Project initiated by Sonia Hasko and maintained by Anirudh Tenneti under the supervision of Tobias Schwaigert.
