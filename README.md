# POMMES v0.1
## Peak-Ordering Measurement for Multi-Elemental Systems

POMMES is a tool written in Python 3 for quantifying order in multi-elemental crystal systems such as perovskite crystals. We aim to deliver a GUI capable of estimating order in multiple crystal sites given ratios between peak intensities in data from X-ray diffraction (XRD).

Currently, we have a Python script into which a user inputs peak ratios in order to obtain an order estimate. Enter intensity ratios between the 001/002, 001/004, and 001/006 peaks into `calc_structure_factor(L1,L2,ratio,color)` and `ratio` will be interpolated in the data for the 00L1/00L2 ratio, plotted in a chosen `color` for an occupancy estimate printed to the console.

Supports EuTa2O6 and SrTa2O6 as of now and will later support other elements.

Project initiated by Sonia Hasko and maintained by Anirudh Tenneti under the supervision of Tobias Schwaigert.
