```@meta
CurrentModule = SeaGap
```

# Static array positioning with sound speed gradients

As discussed in various studies (e.g., [Kido (2007)](https://earth-planets-space.springeropen.com/articles/10.1186/BF03352027), [Honsho et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB017135)), horizontal array displacements has trade-off relationship with horizontal heterogenity of a sound speed structure.
Recent GNSS-A studies has considered a horizontally sloping sound speed structure ([Yokota et al. (2018)](https://link.springer.com/article/10.1007/s11001-018-9362-7), [Yasuda et al. (2017)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL072918), [Honsho et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB017135)).
Here, a static array positioning method considering a horizontally sloping sound speed structure ([Tomita & Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0)) is simply introduced. 

Effect of a horizontally sloping sound speed structure on a travel-time is expressed by two component: shallow and deep gradients (e.g., [Yokota et al., 2018](https://link.springer.com/article/10.1007/s11001-018-9362-7)).

```math
\frac{1}{M\left(\xi_{n,k}\right)} T_{n,k}^{\mathrm{obs}}=\frac{1}{M\left(\xi_{n,k}\right)} {T}^{\mathrm{cal}}\left(\mathbf{p}_k+\color{red}\delta \mathbf{p}\color{black}, \mathbf{u}\left(t_{n},{\bf b}_0\right), v_0\right)+\color{red}C_0\color{black}\left(t_{n}\right) +\color{red}\mathbf{g}_{\mathrm{s}}\color{black} \mathbf{u}^{\mathrm{hor}}\left(t_{n}\right)+\color{red}\mathbf{g}_{\mathrm{d}}\color{black} \mathbf{h}_{n,k}
```
 with

```math
\mathbf{h}_{n,k}=\left(\tan \xi_{n,k} \sin \phi_{n,k}, \tan \xi_{n,k} \cos \phi_{n,k}\right)
```

``\mathbf{g}_{\mathrm{s}}`` is the shallow gradient, and ``\mathbf{g}_{\mathrm{d}}`` is the deep gradient; they individually have EW and NS components.
``\xi`` is a shot angle, and ``\phi`` is an azimuth.
They originally have time flcuatuation, but it is assumed o be constant during a campaign for simplicity.
``C_0`` is NTD corresponding to the average sound speed fluctuation.
The other terms are same with the previous pages ([static array positioning](static.md)).

As long as using single sea-surface platform, the NTD (``C_0``) and the effect of the shallow gradient (``\mathbf{g}_{\mathrm{s}}\mathbf{u}^{\mathrm{hor}}\left(t_{n}\right)``) are strongly trade-off,
This means that the effect of the shallow gradient does not strongly affect horizontal array displacements by the single sea-surface platform case.
Considering this, we simultaneously model the NTD and the effect of the shallow gradient by the superposition of the 3d B-spline functions:

**Observation Equation (0)**
```math
\frac{1}{M\left(\xi_{n,k}\right)} T_{n,k}^{\mathrm{obs}}=\frac{1}{M\left(\xi_{n,k}\right)} {T}^{\mathrm{cal}}\left(\mathbf{p}_k+\color{red}\delta \mathbf{p}\color{black}, \mathbf{u}\left(t_{n},{\bf b}_0\right), v_0\right)+\sum_{j=1}^{J}\color{red}c_j\color{black}\Phi_j(t_n)+\color{red}\mathbf{g}_{\mathrm{d}}\color{black} \mathbf{h}_{n,k}
```

According to the observation equation (0), we can model the sound speed gradient straucture by estimating EW and NS components of the deep gradient (``\mathbf{g}_{\mathrm{d}}``) in addition to the array position and the 3d B-splined NTDs.
You can perform this optimization through Gauss-Newton method by `pos_array_pvg()`.

However, this optimization is quite unstable except acoustic data obtained from spatially well-distributed sea-surface platform (for example, the data of Japan Coast Guard: [Watanabe et al., 2020](https://www.frontiersin.org/articles/10.3389/feart.2020.597532/full)) or obtained at a site with a bunch of seafloor transponders (for example, multi-angled transponders: [Tomita et al., 2019](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-019-1082-y)).

To stabilize this optimization, Honsho et al. (2019) provided a further constaint: the sound speed gradient uniformly exsits from the sea-surface to a certain depth (gradient depth ``D``).
This assumption menas that the deep gradient has a linear relationship with the shallow gradient, and its linear coeffient is the gradient depth as:

```math
{\bf g}_{\rm d}=\frac{D}{2}{\bf g}_{\rm s}
```

Thus, if you somehow obtain the shallow depth, you can estimate the deep gradient by optimizing only the gradient depth:

**Observation Equation (1)**
```math
\frac{1}{M\left(\xi_{n,k}\right)} T_{n,k}^{\mathrm{obs}}=\frac{1}{M\left(\xi_{n,k}\right)} {T}^{\mathrm{cal}}\left(\mathbf{p}_k+\color{red}\delta \mathbf{p}\color{black}, \mathbf{u}\left(t_{n},{\bf b}_0\right), v_0\right)+\sum_{j=1}^{J}\color{red}c_j\color{black}\Phi_j(t_n)+\frac{\color{red}D\color{black}}{2}\mathbf{g}_{\mathrm{s}}\mathbf{h}_{n,k}
```

Meanwhile, because of the trade-off relationship between the temporal fluctuation of NTD and the shallow gradient, it is difficult to grasp the shallow gradient.
Separation of them can be roughly conducted by modeling long-term component of the temporal fluctuation of NTD (e.g., [Honsho et al., 2019](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB017135)).
Here, we model the long-term NTD as the 4th polynomial functions:

**Observation Equation (2)**
```math
\frac{1}{M\left(\xi_{n,k}\right)} T_{n,k}^{\mathrm{obs}}=\frac{1}{M\left(\xi_{n,k}\right)} {T}^{\mathrm{cal}}\left(\mathbf{p}_k+\delta \mathbf{p}, \mathbf{u}\left(t_{n},{\bf b}_0\right), v_0\right)+ \sum_{m=0}^{4} \color{red}\gamma_{m}\color{black}t^m_n + \color{red}\mathbf{g}_{\mathrm{s}}\color{black} \mathbf{u}^{\mathrm{hor}}\left(t_{n}\right) +\frac{D}{2}\color{red}\mathbf{g}_{\mathrm{s}}\color{black}\mathbf{h}_{n,k}
```

Iterative optimization of the observation equations (1) and (2), note that the red terms are the unknown parameters for each equation, can conduct relatively stabler analysis under the sound speed gradient than the observation equation (0) ([Tomita & Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0)).
Moreover, the gradient depth is limited from the sea-surface (zero) to the seafloor.
To flexibly inculude the limitation, we employ a MCMC technique to optimize the observation equations (1) and (2) `pos_array_mcmcpvg()`.
 
To perform `pos_array_mcmcpvg()`, we initially obtain solutions by the static array positioning with horizontally stratified sound speed structure `pos_array_all()`.
Using the solutions of `pos_array_all()` as the initial values, `pos_array_mcmcpvg()` optimize the observation equations (1) and (2) for odd and even MCMC loops, respectively.
Note that we also optimize a scaling factor for each observation equations; the scaling factors express the observational error of each equation and balance the weights of the observation equations (e.g., [Tomita et al., 2021](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JB020991)). 

The details of this MCMC optimization are shown in [Tomita & Kido (2022)](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0).

Optionally, SeaGap includes a developing function `pos_array_mcmcpvgc()` which add some constraints to `pos_array_mcmcpvg()` to obtain solutions more stably.
See Tutorials in detail for this developing function. 

