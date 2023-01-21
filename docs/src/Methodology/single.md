```@meta
CurrentModule = SeaGap
```

# Individual transponder positioning

I order to perform the array positioning, we have to determine initial configuration of a seafloor transponder array (initial positions of the individual seafloor transponders shown in "pxp-ini.xyh") in advance.
Positioning a individual seafloor transponder is performed as similar with the hypocenter determination in Seismology.
Thus, it is important to obtain acoustic shot data with good spatial coverage surrounding a seafloor transponder.
Moreover, to suppress influence of the sound speed fluctuation, it is also important to simultaneously handle multiple seafloor transponder data to estimate precise sound speed fluctuation and to model the sound speed fulcutuation as same with [the static array positioning](static.md).

`pos_single()` is a function to estimate individual seafloor transponder positions and to model the sound speed fluctuation using 3d B-spline functions based on the concept of NTD (Nadir Total Delay).
See [the kinematic array positioning](kinematic.md) page and [Honsho & Kido (2017)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014733) to know NTD in detail.
Then, `pos_single` estimates coefficients for the fixed number of the 3d B-spline functions as well as the ndividual seafloor transponder positions during each campaign through the Gauss-Newton method.

The observation equation for ``n``th acoustic ping transmitted to ``k``th seafloor transponder in `pos_single` is shown as:
```math
\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm obs}_{n,k}=\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm cal}\left({\bf u}(t_n, {\bf b}_0),{\bf p}_k+\color{red}\delta{\bf p}_k\color{black},v_0\right)+\sum_{j=1}^{J}\color{red}c_j\color{black}\Phi_j(t_n)
```

``M`` corresponds to the mapping function, which normalizes a slant travel-time nto a projected travel-time along the perpendicular direction to the seafloor:  
```math
M\left(\xi_{n,k}\right)=\frac{1}{\cos\xi_{n,k}}
```

Note that ``\xi_{n,k}`` is the shot angle from the nadir direction. 

``T^{\rm obs}_{n,k}`` is the observed round-trip travel-time, and ``T^{\rm cal}`` is the calculated travel-time from the assigned parameters: a sea-surface transducer position, a seafloor transponder position, and an underwater sound speed profile.

``{\bf u}(t_n, {\bf b}_0)`` is position of the transducer attached to the sea-surface platform when the shot time is ``t_n``.
The transducer position is transformed from the GNSS antenna position considering the offset between the transducer and the GNSS antenna (``{\bf b}_0`` denoted in "tr-ant.inp") and attitudes of the sea-surface platform.
``{\bf p}_k`` is the initial position of ``k``th transponder, and ``\delta{\bf p}_k``is the unknown paramters of this function, which indicates modification values from the initial positions.
``v_0`` is a sound speed profile (ss\_prof.zv) for calculating travel-times.
``\Phi`` is the 3d B-spline functions using ``J`` bases. NTD was expressed by the 3d B-spline functions with their coefficients of ``c_j``. We estimate the NTD fluctuation by optimizing ``c_j``.

Note that the above ``T^{\rm cal}`` considers difference between the transmitted and recieved sea-surface positions:
```math
T^{\rm cal}\left({\bf u}(t_n),{\bf p}_k,v_0\right)=T^{\rm cal_{\color{blue}{\rm transmit}\color{black}}}\left({\bf u}(t^{\color{blue}{\rm transmit}\color{black}}_n),{\bf p}_k,v_0\right)+T^{\rm cal_{\color{blue}{\rm recieve}\color{black}}}\left({\bf u}(t^{\color{blue}{\rm recieve}\color{black}}_n),{\bf p}_k,v_0\right)
```

