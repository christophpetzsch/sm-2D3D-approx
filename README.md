# Approximate 2D-3D Shape Matching for Interactive Applications

Official repository for the 3DV 2025 paper 'Approximate 2D-3D Shape Matching for Interactive Applications' by Christoph Petzsch, Paul Roetzer, Zorah Lähner and Florian Bernard.


### Installation
This project requires [gptoolbox](https://github.com/alecjacobson/gptoolbox) (click on link and follow instructions on how to install on your machine). You will also need to download [eigen](https://gitlab.com/libeigen/eigen), extract it into `tools/` and rename `eigen-master` into `eigen`. Finally you will be prompted to install further toolboxes by matlab.


### Usage
You can run this code using `experiment_dijkstra.m`. It also contains options to run the code of Lähner et al. and Roetzer et al., as well as options to activate energy normalization and landmark matching. For the latter please use the top bar to enter/exit zoom/rotate/select mode.


### Attribution
When using this code for your own projects please cite the followig:
- Petzsch et al.: Approximate 2D-3D Shape Matching for Interactive Applications, 3DV 2025
- Roetzer et al.: Conjugate Product Graphs for Globally Optimal 2D-3D Shape Matching, CVPR 2023
- Lähner et al.: Efficient Globally Optimal 2D-to-3D Deformable Shape Matching, CVPR 2016


### License
This repo is licensed under MIT license.


