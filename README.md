# k-wSST-PDA-WindFarm

<!--  [![Compatibility: OFver](https://img.shields.io/badge/Compatible_with-OpenFOAM.v2112-lightblue.svg)]()  -->
[![Paper: RENE](https://img.shields.io/badge/Reference-Paper-red.svg)](https://doi.org/10.48550/arXiv.2405.04906)
[![Paper: Author](https://img.shields.io/badge/Author-green.svg)](https://sites.google.com/view/zehtabiyan/home)

<!-- [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)   -->

# A progressive data-augmented RANS model for enhanced wind-farm simulations 
A progressive data-augmented version of $k-\omega\text{SST}$ model, called $k-\omega\text{SST-PDA}$, integrated with OpenFOAM for simulation of turbine wakes and power losses in wind farms, developed by [Amarloo and Zehtabiyan-Rezaie and Abkar (2024)](https://doi.org/10.48550/arXiv.2405.04906) at the Fluid Mechanics and Turbulence research group at Aarhus University, Denmark. 

# Description
Here, we propose a progressive data-augmentation approach to enhance the accuracy of Reynolds-averaged Navier-Stokes models to represent flow characteristics within wind-turbine and wind-farm wakes. We first incorporate the turbine-induced forces in the turbulent kinetic energy equation of the widely used $k−\omega\text{SST}$ model as proposed in this [work](https://doi.org/10.1016/j.renene.2021.08.012). Afterward, we utilize data from large-eddy simulations to progressively enhance the Reynolds-stress prediction of this baseline model, accurately capturing the evolution of eddy viscosity in the wake, as well as the emergence of secondary flows. To read more about the second kind of secondary flows in the context of wind-turbine wake, see [our paper](https://doi.org/10.1063/5.0203068). 

# Target platform
The code has been rigorously tested and verified to be fully compatible with OpenFOAM v-2312, ensuring its smooth integration and reliable performance with this specific release.

<!-- # Author
[Navid Zehtabiyan-Rezaie](https://sites.google.com/view/zehtabiyan/home)-->

# How to set the model
1- Download the source code.
 
2- Copy the folder _OF_kOmegaSSTWindFarm_ to your library directory, and compile the turbulence model by the following command: 

`wmake`

3- Copy the folder _testCase_ to your run directory and execute.

4- To learn more about the momentum and turbulent kinetic energy sources linked to the turbine forces, see this repository: https://github.com/AUfluids/k-epsilon-Sk.

# Evaluation of $k-\omega\text{SST-PDA}$ model's performance
In our latest [publication](https://doi.org/10.48550/arXiv.2405.04906), we apply the optimized $k−\omega\text{SST-PDA}$ model to several wind-farm cases with distinct layouts and conduct a comparative analysis focusing on the obtained quantities such as normalized streamwise velocity deficit, turbulence intensity, and power output. We examine the success rate of the augmented model in predicting the secondary flows in the wake region. Our comparisons and validations demonstrate the superior performance of the progressive data-augmented model over the standard version in all cases considered.

<img src="https://github.com/AUfluids/k-wSST-PDA-WindFarm/blob/main/testCase/Performance.png" width="900" height="450" alt="Normalized power of turbines in three validation cases">

  
# How to cite
Please, cite this library as:
```
@misc{amarloo2024PDA,
      title={A progressive data-augmented RANS model for enhanced wind-farm simulations}, 
      author={Ali Amarloo and Navid Zehtabiyan-Rezaie and Mahdi Abkar},
      year={2024},
      eprint={2405.04906},
      archivePrefix={arXiv},
      primaryClass={physics.flu-dyn},
      url={https://arxiv.org/abs/2405.04906}, 
}
```
