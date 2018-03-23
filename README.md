# Listeria_spread_simulations
Monte Carlo simulations of Listeria monocytogenes cell-cell spread

All the simulations plot (1) the number of bacteria, (2) the second moment, (3) the area of the convex hull around the data,
and (4) the radial distance to the center of mass, as functions of time. The growth rate (k) can be calculated by fitting
the number of bacterial to an exponential function. The diffusion coefficient (D) can be calculated from the slope of the
second moment given that it equals 4Dt for a two-dimensional diffusive process. The speed of the bacterial front
of a reaction-diffusion equation (diffusion + growth) equals the sqrt(4Dk) as the limit approaches infinity. This number
is approximately similar to the slope of the radial distance to the center of mass.

**Note:** you can uncomment lines 7-19 of **simConvexHullNEW.m** to watch the progression of the simulation. It will be plotted
according to the variable called **plotevery** in line 19 of the main program.

**Gaussian_sim.m:** simulation of Listeria spread where the size of the hop is sampled from a Normal distribution
with mean = 0 and std = sqrt(2Dt). The file **All bacteria + time.mov** is a sample movie created by this simulation. If
two diffusion coefficients are used simultaneously, then a movie such as **TwoDMovie + Time.avi** is produced.

**HeavyTail_sim.m:** simulation of Listeria spread where the size of the hop is sampled (using inverse transform sampling)
from a custom heavy-tailedp polynomial 1/u^alpha, where alpha determines the heaviness of the distribution's tail. The
distribution's second moment converges only when alpha > 3.

**Example_Foci.pdf** shows an example of the output of HeavyTail_sim.m (left) and Gaussian_sim.m (right).

**GaussianExtrusionTWO.m:** simulation of Listeria spread where the bacteria move Gaussianly, but the process is governed
by two diffusion coefficients. The ratio of the diffusion coefficients is set by **cFactor**. In addition, the variable
**hostCell** mimics the extrusion of epithelial cells from the intestinal villi. In the simulation, host cell extrusion
events lead to the death of bacteria inside the radius of the host cell. The simulation keeps track of the number of
extruded bacteria as well as the number of bacteria spreading through the epithelium.
