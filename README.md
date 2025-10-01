# Geometric Spatio-Spectral Total Variation for Hyperspectral Image Denoising and Destriping

This is a demo code of the proposed method in the following reference:

S. Takemoto and S. Ono,
``Geometric Spatio-Spectral Total Variation for Hyperspectral Image Denoising and Destriping.''

Update history:
1, Oct. 2025: v1.0 

For more information, see the following

- Preprint paper: 

## How to use
1. **Downloding datasets
 - Jasper Ridge: https://rslab.ut.ac.ir/data
 - Pavia University: https://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes#Pavia_Centre_and_University
2. **Setting parameters**
 - Choose the image (JasperRidge or PaviaUniversity)
 - Adjust the parameters
   - `params.rho`: parameter for the radii of the noise terms
   - `params.omega`: balancing parameter between first- and second-order differences
   - `params.stopcri`: Stopping criterion
   - `params.maxiter`: Maximum number of iterations
   - `params.disprate`: Period to display intermediate results
 - Set as `use_GPU` = 1 if you use GPU.

3. Run ```main_GeoSSTV.m```


## Our Reference
If you use this code, please cite the following paper:

```

```

