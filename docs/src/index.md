```@meta
CurrentModule = SeaGap
```

# SeaGap

SeaGap is written in Julia for fast calculation. SeaGap is designed to perform "text-based processing" that SeaGap reads data from the input text file, performs calculation, and then writes the output text file instead of holding the calculation results as variables in Julia.

GNSS-Acoustic (GNSS-A) technique is a combination technique of GNSS positioninig on a sea-surface platform and acoustic ranging between a sea-surface transducer for seafloor geodesy, which is contrived by Dr. Spiess from the Scripps Institution of Oceanography in 1980s. The GNSS-A technique has been developed by various researchers. Altanative tool for GNSS-A positioning is [GARPOS](https://github.com/s-watanabe-jhod/garpos) develped by [Watanabe et al. (2020)](https://www.frontiersin.org/articles/10.3389/feart.2020.597532/full).

SeaGap provides various functions to easily perform GNSS-A positioning, post-processing for the positioning results, and visualization of the results. Moreover, SeaGap covers various types of positioning: classic kinematic array positioning (e.g., [Kido et al. 2006](https://earth-planets-space.springeropen.com/articles/10.1186/BF03351996)), static array positioning (e.g., [Tomita and Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0)), MCMC-base static array positioning (e.g., [Tomita and Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0)), static positioning for an individual transponder, and static array positioning with estimation of an offset between a sea-surface GNSS antenna and a sea-surfacetransducer.

You should refer  [Tomita and Kido, 2022](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01740-0) for the citation of this software (note that an article  for this software is preparing). 


**Contents:**
```@contents
Pages = [
      "Dataformat/intro.md",
      "Methodology/intro.md",
      "Tutorials/intro.md",
      "Others/intro.md",
]
Depth = 2
```

