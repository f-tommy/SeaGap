```@meta
CurrentModule = SeaGap
```

# Overview

In this page, overview of the GNSS-Acoustic (GNSS-A) positioning methods is introduced.

GNSS-A positioning was contrived by [Spiess (1985)](https://ieeexplore.ieee.org/document/4072330), and the geodetical results were firstly given by [Spiess et al. (1998)](https://www.sciencedirect.com/science/article/abs/pii/S0031920198000892).
Then, it has been developed by various researchers, especially in Japan (Japan Coast Guard, Nagoya Univ., and Tohoku Univ.).

Each GNSS-A site is composed of multiple seafloor transponders (minimum: 3), which is called a seafloor transponder array.
In the GNSS-A positioning, we first determine initial positions of the seafloor transponders (corresponding to "pxp-ini.xyh" in SeaGap).
Then, we estimate an array displacement of the seafloor transponders from the initial positions assuming the seafloor transponder array rigidly moves among the repeated GNSS-A campaigns.
This estimation is performed by estimating a common position change among the seafloor transponders for each campaign.

The most accurate technique to determin the initial positions of the seafloor transponders is the estimation using various campaign data (e.g., [Honsho & Kido, 2017](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JB014733)).
The SeaGap do not have such an elaborated function, but has a function to determine positions of the individual seafloor transponders using single campaign data [Individual transponder positioning `pos_single()`](single.md).
It is a simple but powerful way that we first estimate positions of the individual seafloor transponders from single campaign data and then average them from multiple campaigns as the initial positions.

The array positioning is separated into two approaches: kinematic and static positioning.
[The kinematic positioning approach](kinematic.md) estimates an array displacement using very short-term positioning data, and [the static positioning apprach](static.md) estimates an array displacement using a whole campiagn data.
The charactaristics of each approach are shown in their individual pages.



