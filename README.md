# [NeurIPS 2023 (Spotlight)] “Robust Model Reasoning and Fitting via Dual Sparsity Pursuit” 


### Requirements
Successfully test with Matlab 2016b! 

VLFeat: https://www.vlfeat.org/ （if vlfeat is not workable, please redownload then decompress it into subFun/）

# Run the Code

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

# Results
### Unknown model reasoning results:

qualitative result   |  quantitative result
:-------------------:|:--------------------------------------:
![](https://github.com/StaRainJ/DSP/blob/main/fig/FigMatchResults.png)  |  ![](https://github.com/StaRainJ/DSP/blob/main/fig/TabModelReasoning.png)

### Exact model fitting results:
[hard datasets]

![image](https://github.com/StaRainJ/DSP/blob/main/fig/FigDatasets1.png)

[easy datasets]

![image](https://github.com/StaRainJ/DSP/blob/main/fig/FigDatasets2.png)

# Citation
### Please kindly cite this paper if you use the code in this repository as part of a published research project.

 Jiang Xingyu, Ma Jiayi. Robust Model Reasoning and Fitting via Dual Sparsity Pursuit. NeurIPS 2023 
```bash
 @inproceedings{jiang2023dsp,
	title={Robust Model Reasoning and Fitting via Dual Sparsity Pursuit},
	author={Xingyu Jiang, Jiayi  Ma},
        booktitle={Advances in Neural Information Processing Systems},
	year={2023}
}
```
