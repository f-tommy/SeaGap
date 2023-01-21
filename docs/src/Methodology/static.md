```@meta
CurrentModule = SeaGap
```

# Static array positioning

## Overview

Static array positioning estimates an array displacement using a whole campaign data for each site.
This technique has been developed by Japan Coast Guard (e.g., [Fujita et al., 2006](https://earth-planets-space.springeropen.com/articles/10.1186/BF03351923); [Sato et al., 2013](https://link.springer.com/article/10.1007/s00190-013-0649-9); [Watanabe et al., 2020](https://www.frontiersin.org/articles/10.3389/feart.2020.597532/full)), Nagoya Univ. (e.g., [Ikuta et al., 2008](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2006JB004875); [Kinugasa et al., 2020](https://progearthplanetsci.springeropen.com/articles/10.1186/s40645-020-00331-5)), and recently Tohoku Univ. ([Honsho & Kido, 2017](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014733); [Tomita & Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0)). 

In the static arrat positioning approach, a temporal sound speed fluctuation is generally modeled by a kind of flexible functions, such as 3d B-spline functions (e.g., [Honsho & Kido, 2017](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014733)). 
The static array positioning in SeaGap `pos_array_all` modeled a temporal sound speed fluctuation as the NTD (Nadir Total Delay; see [the kinematic array positioning](kinematic.md) page for detail) temporal variation and expresses it by 3d B-spline functions. 
Then, `pos_array_all` estimates coefficients for the fixed number of the 3d B-spline functions as well as the array position averagely during each campaign through the Gauss-Newton method.

## Basic static array positioning method

The observation equation for ``n``th acoustic ping transmitted to ``k``th seafloor transponder in `pos_array_all` is shown as:
```math
\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm obs}_{n,k}=\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm cal}\left({\bf u}(t_n, {\bf b}_0),{\bf p}_k+\color{red}\delta{\bf p}\color{black},v_0\right)+\sum_{j=1}^{J}\color{red}c_j\color{black}\Phi_j(t_n)
```

``M`` corresponds to the mapping function, which normalizes a slant travel-time nto a projected travel-time along the perpendicular direction to the seafloor:  
```math
M\left(\xi_{n,k}\right)=\frac{1}{\cos\xi_{n,k}}
```

Note that ``\xi_{n,k}`` is the shot angle from the nadir direction. 

``T^{\rm obs}_{n,k}`` is the observed round-trip travel-time, and ``T^{\rm cal}`` is the calculated travel-time from the assigned parameters: a sea-surface transducer position, a seafloor transponder position, and an underwater sound speed profile.

``{\bf u}(t_n, {\bf b}_0)`` is position of the transducer attached to the sea-surface platform when the shot time is ``t_n``.
The transducer position is transformed from the GNSS antenna position considering the offset between the transducer and the GNSS antenna (``{\bf b}_0`` denoted in "tr-ant.inp") and attitudes of the sea-surface platform.
``{\bf p}_k`` is the position of ``k``th transponder.
``\delta{\bf p}`` is the unknown parameter corresponding to the array displamcent, and the 3 components (EW, NS, UD) of the array position are estimated.
``v_0`` is a sound speed profile (ss\_prof.zv) for calculating travel-times.
``\Phi`` is the 3d B-spline functions using ``J`` bases. NTD was expressed by the 3d B-spline functions with their coefficients of ``c_j``. We estimate the NTD fluctuation by optimizing ``c_j``.

Note that the above ``T^{\rm cal}`` considers difference between the transmitted and recieved sea-surface positions:
```math
T^{\rm cal}\left({\bf u}(t_n),{\bf p}_k+\delta{\bf p},v_0\right)=T^{\rm cal_{\color{blue}{\rm transmit}\color{black}}}\left({\bf u}(t^{\color{blue}{\rm transmit}\color{black}}_n),{\bf p}_k+\delta{\bf p},v_0\right)+T^{\rm cal_{\color{blue}{\rm recieve}\color{black}}}\left({\bf u}(t^{\color{blue}{\rm recieve}\color{black}}_n),{\bf p}_k+\delta{\bf p},v_0\right)
```


## Optimization of number of the basis functions

In SeaGap, the bases of the 3d B-spline function are deployed with temporally constant intervals, and the number of them can be optionally provided. 
The number of the bases of the 3d B-spline function directly affects temporal smoothness of the modeled NTD fluctuation so that it should be optimized depending on the observational data.

Various criteria are proposed to optimize selection of unknown paramters.
SeaGap supports to calculate two major criterion "AIC (Akaike's information criterion)" and "BIC (Bayesian Information Criterion)" when performing the above basic static array positioning.
The both criterion have their own charactersitics.
The BIC is a model for situations where the error range of the maximum likelihood estimate is extremely small compared to the value of the parameter, which means that significant parameters are easily distinguished from non-significant ones; in contrast, AIC focuses on the treatment of parameters that are only marginally significant, and pursues modeling possibilities to the point where they are almost buried in error effects ([Akaike, 1996](https://orsj.org/wp-content/or-archives50/pdf/bul/Vol.41_07_375.pdf) [in Japanese]).
As for the case that we optimize the number of the bases of the 3d B-spline function, AIC tends to support large number because AIC is sensitive to a local temporal trend which is almost buried in error effects.
Thus, we genrally employ BIC for this optimization ([Tomita & Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0))

In SeaGap, `pos_array_all()` function automatically returns the AIC and BIC values for a given number of the bases, and `pos_array_all_AICBIC()` function returns the AIC and BIC values for various number of the bases.

## Optimization of an offset between a transducer and a GNSS antenna

In SeaGap, `pos_array_TR()` function optimizing an offset between a transducer and a GNSS antenna is prepared as a kind of the static array positioning approach.
The optimization of the offset is essential when it cannot be accurately measured by the other technique (such as optical ranging).
Importance of this procedure is documented in [Honsho et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB017135).

The observation equation for ``n``th acoustic ping transmitted to ``k``th seafloor transponder in `pos_array_TR` is shown as:
```math
\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm obs}_{n,k}=\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm cal}\left({\bf u}(t_n, {\bf b}_0+\color{red}\delta {\bf b}\color{black}),{\bf p}_k+\color{red}\delta{\bf p}\color{black},v_0\right)+\sum_{j=1}^{J}\color{red}c_j\color{black}\Phi_j(t_n)
```

``\delta {\bf b}`` is the modification value from the initial offset written in "tr-ant.inp".
This value is obtained from each campaign data at single site; thus, averaging the modification values for all campaign data using the same configuration of the GNSS antenna and the transducer is required to obtain accurate value.


