#  Code of “Robust Model Reasoning and Fitting via Dual Sparsity Pursuit” NeurIPS 2023 (Spotlight)
### Requirements
Successfully test on Matlab 2016b! 

VLFeat: https://www.vlfeat.org/

# Run the code

### for exact model fitting via DSP
Test for two-view geometric model (H or F):
```bash
run Demo4KnownModel.m
```

### for unknown model reasoning
Test for 2D point data:
```bash
run Demo4point.m
```
Test for two-view data:
```bash
run Demo4UnknownModel.m
```

# Method Description
### Illustration of sparse subspace recovery:

![image](https://github.com/StaRainJ/DSP/blob/main/fig/Fig1.png)

### Data Embedding of Existing Methods and our DSP:

![image](https://github.com/StaRainJ/DSP/blob/main/fig/TabDataEmbedding.png)


# Please kindly cite this paper if you use the code in this repository as part of a published research project.

 Jiang Xingyu, Ma Jiayi. Robust Model Reasoning and Fitting via Dual Sparsity Pursuit. NeurIPS 2023 
