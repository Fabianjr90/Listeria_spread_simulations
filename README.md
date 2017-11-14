# Listeria_spread_simulations
Monte Carlo simulations of Listeria monocytogenes cell-cell spread

All the simulations plot (1) the number of bacteria, (2) the second moment, (3) the area of the convex hull around the data,
and (4) the radial distance to the center of mass, as functions of time. The growth rate (k) can be calculated by fitting
the number of bacterial to an exponential function. The diffusion coefficient (D) can be calculated from the slope of the
second moment given that it equals 4Dt for a two-dimensional diffusive process. The speed of the bacterial front
of a reaction-diffusion equation (diffusion + growth) equals the sqrt(4Dk) as the limit approaches infinity. This number
is approximately similar to the slope of the radial distance to the center of mass.

**Gaussian_sim.m:** simulation of Listeria spread where the size of the hop is sampled from a Normal distribution
with mean = 0 and std = sqrt(2Dt).

**HeavyTail_sim.m:** simulation of Listeria spread where the size of the hop is sampled (using inverse transform sampling)
from a custom heavy-tailedp polynomial 1/u^alpha, where alpha determines the heaviness of the distribution's tail. The
distribution's second moment converges only when alpha > 3.
