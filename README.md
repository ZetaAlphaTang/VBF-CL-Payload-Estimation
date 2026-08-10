# VBF-CL Experimental Code

This repository contains the MATLAB code and experimental data used for the three indoor UAV suspended-payload experiments reported in the paper.

## Folder Structure

- `Experiment_1/`: hovering experiment
- `Experiment_2/`: straight-line motion experiment
- `Experiment_3/`: figure-eight maneuver experiment

Each folder contains the corresponding `.mat` data files, the main script `data_processing.m`, and the filtering functions used in the comparison.

## Requirements

- MATLAB
- No additional toolbox is required for the filtering code. MATLAB video export is used to generate the optional trajectory animation.

## How to Run

Open MATLAB, set the current folder to one experiment directory, and run:

```matlab
data_processing
```

For example:

```matlab
cd('Experiment_1')
data_processing
```

Repeat the same procedure for `Experiment_2` and `Experiment_3`.

## Outputs

The script prints the RMSE statistics in the MATLAB Command Window and generates comparison figures for payload position and velocity estimation. It also creates a 3-D trajectory animation file named:

- `Experiment1_3D_Trajectory_Animation.mp4`
- `Experiment2_3D_Trajectory_Animation.mp4`
- `Experiment3_3D_Trajectory_Animation.mp4`

depending on the selected experiment.

## Implemented Filters

The comparison includes EKF, EKF without input, PL-KF, PL-KF without input, SHF, SHF without input, VBF without input, and the proposed VBF-CL method.

## Notes

The scripts use relative paths, so each `data_processing.m` file should be executed from its own experiment folder. If video generation is not needed, the call to `create_figure1_trajectory_video` in `data_processing.m` can be commented out.

## Citation

If you use this code, please cite the associated paper:

J. Tang, Y. Jin, Y. Li, Y. Xu, X. Yang, and B. Zhu, "Cable-Length-Aided Variational Bayesian Filtering for Suspended-Payload Motion Estimation Using Low-Cost Onboard Sensors."
