This is a MatLab/C++ library of [1].

The code supports arbitrary number of regions and inclusion and exclusion constraints.
The syntax is described in "example.m".

Requirements 
---
1. MATLAB and a c++ compiler.

Getting started
---
1. Run mex --setup.
2. Run example.m in MATLAB, the mex'ed c++ code compiles automatically.

All functionally is described in example.m

Third party software
---
* To calculate the correct weights for the regularization the code uses
  Sphere voronoi by John Burkardt [2].
  

* If you use the solver based on lagrangian duality, the code uses the max/flow min solver
	By Yuri Boykov and Vladimir Kolmogorov [3] which reuses flow as described [4].

* If you use roof duality, the code uses the QPBO software [5]

More resources
---
Martin Rykfors has written a more efficient solver 
if you problem only have inclusions and you are satisfied
with standard connectivities.
[Code](https://github.com/MartinRykfors/MultiRegion)
More details in his [Master thesis](http://www.maths.lth.se/vision/education/pages/Rykfors12/exjobb.pdf)

References
----------

1. [An Efficient Optimization Framework for Multi-Region Segmentation based on Lagrangian Duality](http://www.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=531)
	 Johannes Ulén, Petter Strandmark and Fredrik Kahl. 

2. [Sphere voronoi](http://people.sc.fsu.edu/~jburkardt/m_src/sphere_voronoi/sphere_voronoi.html)
John Burkardt

3 [An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision]
	Yuri Boykov and Vladimir Kolmogorov.

4. [Dynamic graph cuts for efficient inference in markov random fields.]
	Kohli, Pushmeet, and Philip HS Torr.

5.[Optimizing binary MRFs via extended roof duality]
	C. Rother, V. Kolmogorov, V. Lempitsky, and M. Szummer.