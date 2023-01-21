```@meta
CurrentModule = SeaGap
```

# Kinematic array positioning


Kinematic array positioning is the classical GNSS-Acoustic positioning techqniue, which firstly contrived by [Spiess et al. (1998)](https://www.sciencedirect.com/science/article/abs/pii/S0031920198000892). In the classical Scrips transponder system and the Tohoku Univ.'s transponder system, all seafloor transponders reply to a sea-surface acoustic ping; while a sea-surface transducer calls each seafloor transponder in the Japan Coast Guard's system. The former system can collect replies from all transponders simultaneously; thus, kinematic array positioning, which estimates a horizontal array position for each acoustic ping, has been developed by the Scrips institute and Tohoku University.

[Kido et al. (2006)](https://earth-planets-space.springeropen.com/articles/10.1186/BF03351996) and [Kido et al. (2008)](https://earth-planets-space.springeropen.com/articles/10.1186/BF03352785) improved the kinematic array positioning technique by introducing the concept of NTD (Nadir Total Delay) that is analogy of ZTD (Zenith Total Delay) in the GNSS positioning technique. This concept enable us to accurately model a temporal sound speed fluctuation. The detailed explanation of NTD is written in [Honsho & Kido (2017)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014733) and [Tomita et al. (2019)](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-019-1082-y).

Recently, the kinematic array positioning has been further developed by [Tomita et al. (2019)](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-019-1082-y); an extended Kalman filter was employed to the kinematic array positioning. Moreover, this paper also found that kinematic array positioning in vertical component can be performed at a site with multi-angled transponder deployment. But, the kinematic array positioning using the extended Kalman filter is not equipped in the current version of SeaGap.

In SeaGap, the kinematic positioning function (`pos_array_each`) based on [Kido et al. (2006)](https://earth-planets-space.springeropen.com/articles/10.1186/BF03351996) and [Kido et al. (2008)](https://earth-planets-space.springeropen.com/articles/10.1186/BF03352785) is provided. A kinematic array position can be estimated for each shot group by Gauss-Newton method.

The observation equation for ``n``th acoustic ping transmitted to ``k``th seafloor transponder in `pos_array_each` is shown as:
```math
\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm obs}_{n,k}=\frac{1}{M\left(\xi_{n,k}\right)}T^{\rm cal}\left({\bf u}(t_n,{\bf b}_0),{\bf p}_k+\color{red}\delta{\bf p_n}\color{black},v_0\right)+\color{red}C\color{black}
```

``M`` corresponds to the mapping function, which normalizes a slant travel-time nto a projected travel-time along the perpendicular direction to the seafloor:
```math
M\left(\xi_{n,k}\right)=\frac{1}{\cos\xi_{n,k}}
```

Note that ``\xi_{n,k}`` is the shot angle from the nadir direction. 

``T^{\rm obs}_{n,k}`` is the observed round-trip travel-time, and ``T^{\rm cal}`` is the calculated travel-time from the assigned parameters: a sea-surface transducer position, a seafloor transponder position, and an underwater sound speed profile.

``{\bf u}(t_n, {\bf b}_0)`` is position of the transducer attached to the sea-surface platform when the shot time is ``t_n``.
The transducer position is transformed from the GNSS antenna position considering the offset between the transducer and the GNSS antenna (``{\bf b}_0`` denoted in "tr-ant.inp") and attitudes of the sea-s
urface platform.
``{\bf p}_k`` is the position of ``k``th transponder. ``\delta{\bf p_n}`` is the unknown parameter corresponding to the array displamcent, and this value is estimated for each shot group by `pos_array_each`. Although the position of a transponder is expressed in 3 components, the array displcement is estimated only in the horizontal components; threfore, ``\delta{\bf p}=\left(\delta p_{\rm EW},\delta p_{\rm NS},0\right)``. ``v_0`` is a sound speed profile (ss\_prof.zv) for calculating travel-times. ``C`` is NTD expressing temporal perturbation of the sound speed, which is also the unknow paramter as well as the array displacement.

Note that the above ``T^{\rm cal}`` considers difference between the transmitted and recieved sea-surface positions:
```math
T^{\rm cal}\left({\bf u}(t_n),{\bf p}_k+\delta{\bf p},v_0\right)=T^{\rm cal_{\color{blue}{\rm transmit}\color{black}}}\left({\bf u}(t^{\color{blue}{\rm transmit}\color{black}}_n),{\bf p}_k+\delta{\bf p},v_0\right)+T^{\rm cal_{\color{blue}{\rm recieve}\color{black}}}\left({\bf u}(t^{\color{blue}{\rm recieve}\color{black}}_n),{\bf p}_k+\delta{\bf p},v_0\right)
```
The expression of ``{\bf b}_0`` is omitted in the above equation.

In `pos_array_each`, you can assign the shot group number to each ping in the input file (obsdata.inp). Simultaneously using travel-times for each shot group, `pos_array_each` estimates the horizontal array displacement (``\delta{\bf p}``) and NTD (``C``).

