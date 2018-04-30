# Bayesian_Optimization_Material_design

## Accelerated Design of Layered Materials with Bayesian Optimization

#### 1 ```predict_maxval.py:```
*A Gaussian Process Regression model to predict band gap.* 

#### 2 ```predict_structure.py:```
*A Gaussian Process Regression model to predict conduction band maxima, valance band minima and Thermoelectic EFF function.*

#### 3 ```Bayesian_opt.py:``` 
*A Bayesian Optimization model to discover material with optimal property witn minimum structure evalulation*

*Dataset* <br /> 
```N_doped_EFF_max.txt:``` EFF(Electronic fitness function) maximum value for n-doped 3-layer heretero-structure <br /> 
```P_doped_EFF_max.txt:``` EFF(Electronic fitness function) maximum value for p-doped 3-layer heretero-structure  <br /> 
```3-layer-band_gap.txt:``` Maximum Band gap for 3-layer heretero-structure <br /> 

*Input Paramaters:* <br /> 
Inside the code Bayesian_opt.py, we have <br /> 
Nruns = 1                                *```total number of Bayesian  Optimization runs```* <br /> 
train_test_split=0.10                    *```initial sampled data in a given Bayesian  Optimization run```* <br /> 

*To find n-doped 3-ayer hetero-structure with optimal EFF value* <br /> 
Run ```python3.6 Bayesian_opt.py N_doped_EFF_max ``` <br /> 
*To find p-doped 3-ayer hetero-structure with optimal EFF value* <br /> 
Run ```python3.6 Bayesian_opt.py P_doped_EFF_max ``` <br /> 
*To find maximum band gap* <br />
Run ```python3.6 Bayesian_opt.py 3-layer-band_gap.txt``` <br /> 










