# Bayesian_Optimization_Material_design

## Accelerated Design of Layered Materials with Bayesian Optimization

#### 1 ```predict_maxval.py:```
*A Gaussian Process Regression model to predict band gap.* 

predict_maxval.py: creates a gaussian process regression model using x% of data from 3-layer-band_gap.txt. After buiding the model predcits the band gap of the remaining (1-x%) data as test set. <br />

*Dataset* <br /> 
```3-layer-band_gap.txt``` CBM and VBM 3-layer heretero-structure <br /> 
```xaxisvalue.txt```   Contains the value of the wave vector 

*Input Paramaters:* <br /> 
inputfile_name="3-layer-band_gap.txt"    * ###file name of the input data * <br />
train_test_split=0.60                    * ###split between training and test set* <br />
Nrun = 1    <br />

To run the program: ```predict_maxval.py ``` <br /> 

#### 2 ```predict_structure.py:```
*A Gaussian Process Regression model to predict conduction band maxima, valance band minima and Thermoelectic EFF function.*

predict_structure.py build a gaussian process (GP) regression model for condunction band minima and valance band maximum. It takes 3-layer-band_structure.txt as inout data and split in into training and test set. GP model is build using training set and the images of prediced condunction band minima and valance band maximum of the test set is saved in folder Bandstructure.

*Dataset* <br /> 
```"3-layer-band_structure.txt``` Band gap for 3-layer heretero-structure <br /> 

*Input Paramaters:* <br /> 
inputfile_name="3-layer-band_structure.txt"    * ###file name of the input data * <br />
train_test_split=0.60                    * ###split between training and test set* <br />
Nrun = 1    <br />

To run the program: ```predict_structure.py ``` <br /> 

#### 3 ```Bayesian_opt.py:``` 
*A Bayesian Optimization model to discover material with optimal property witn minimum structure evalulation*

*Dataset* <br /> 
```N_doped_EFF_max.txt ```  EFF (Electronic fitness function) maximum value for n-doped 3-layer heretero-structure <br /> 
```P_doped_EFF_max.txt ```  EFF (Electronic fitness function) maximum value for p-doped 3-layer heretero-structure  <br /> 
```3-layer-band_gap.txt ``` Maximum Band gap for 3-layer heretero-structure <br /> 

*Input Paramaters:* <br /> 
Inside the code Bayesian_opt.py, we have <br /> 
Nruns = 1                  * ### total number of Bayesian  Optimization runs* <br /> 
train_test_split=0.10      * ### initial sampled data in a given Bayesian  Optimization run* <br /> 

*To find n-doped 3-ayer hetero-structure with optimal EFF value* <br /> 
Run ```python3.6 Bayesian_opt.py N_doped_EFF_max ``` <br /> 
*To find p-doped 3-ayer hetero-structure with optimal EFF value* <br /> 
Run ```python3.6 Bayesian_opt.py P_doped_EFF_max ``` <br /> 
*To find maximum band gap* <br />
Run ```python3.6 Bayesian_opt.py 3-layer-band_gap.txt``` <br /> 










