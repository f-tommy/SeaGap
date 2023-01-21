```@meta
CurrentModule = SeaGap
```

# Travel-time calculation

## Overview

In SeaGap, travel-time calculation is performed by approximate calculation contrived by [Tomita and Kido (2022)](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0). The approximate travel-time can be speedly obtained, and its accuracy is comparable to that of a general exact travel-time.    
The general features of both the exact and approximate travel-time calculation are shown below; however, as for details of the travel-time calculation, please refer the above paper.

## Exact travel-time calculation

Assuming a horizontally stratified sound speed structure, the shooting method is performed to calculate travel-times. In the calculation, Snell’s law in a spherical frame is considered. Using a local radius depending on latitude (Gaussian radius: [Chadwell and Sweeney 2010](https://www.tandfonline.com/doi/abs/10.1080/01490419.2010.492283?journalCode=umgd20), travel-times can simply be calculated.

The sound speed structure is expressed as a horizontally-layred structure, and a ray parameter is preserved among the layers based on the Snell's law. Optimizing the ray parameter by non-linear least squared method to minimize the supposed and calculated angular distance between a sea-surface transmitting point and a seafloor recieving point, we can obtain an accurate travel-time (``T_{\rm Exact}``).

## Approximate travel-time calculation

To calculate an approximate travel-time (``T_{\rm Appr}``}), a travel-time in a sphrical Earth without considering the acoustic-line bending (``T_{\rm Sphere}``). This can be easily calculated. Then, we add correction terms ``f`` to ``T_{\rm Sphere}``; this means ``T_{\rm Exact}{\approx}T_{\rm Appr}=T_{\rm Sphere}+f``. The correction terms are expressed by 8th order polynomial functions, and these terms are optimized in advance.



