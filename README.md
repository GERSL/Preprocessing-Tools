# Preprocessing-Tools
Preprocessing tools for Landsat data: BRDF (Bidirectional Reflectance Distribution Function) and Topographic Corrections.

Please contact Shi Qiu (shi.qiu@uconn.edu), Rong Shang (shangrong@uconn.edu), and Zhe Zhu (zhe@uconn.edu) at Department of Natural Resources and the Environment, University of Connecticut if you have any questions.

------------

The BRDF correction is to use the c-factor approach [(Roy, D. P. et al., 2016)](https://doi.org/10.1016/j.rse.2016.01.023) based on the RossThick-LiSparse-R BRDF model [(Schaaf, Crystal B., et al. 2002)](https://doi.org/10.1016/S0034-4257(02)00091-3).

The SCS correction is equivalent to projecting the sunlit canopy from the sloped surface to the horizontal surface in the direction of illumination [(Gu, D. et al., 1998)](https://doi.org/10.1016/S0034-4257(97)00177-6).

The SCS+C model is based on the same SCS model, but it integrates a semi-empirical parameter (C) that can significantly reduce the overcorrection caused by the scattered radiation from the source of illumination [(Soenen, S. A. et al., 2005)](https://10.1109/TGRS.2005.852480). Also, the computing efficiency was improved by a  stratified sampling approach [(Qiu, S., et al., 2017)](https://doi.org/10.1016/j.rse.2017.07.002).

The IC model was proposed to remove the dependency of the reflectance on the cosine of the local solar incidence angle (cosi) based on the same linear regression shown [(Tan, B. et al., 2013)](https://doi.org/10.1016/j.rse.2013.05.013).

-------------
If using those functions, please cite the following papers:

paper 1: Qiu, S., Lin, Y., Shang, R., Zhang, J., Ma, L. and Zhu, Z., 2019. Making Landsat Time Series Consistent: Evaluating and Improving Landsat Analysis Ready Data. Remote Sensing, 11(1), p.51.[https://doi.org/10.3390/rs11010051](https://doi.org/10.3390/rs11010051).

paper 2: Qiu, S., He, B., Zhu, Z., Liao, Z. and Quan, X., 2017. Improving Fmask cloud and cloud shadow detection in mountainous area for Landsats 4–8 images. Remote Sensing of Environment, 199, pp.107-119. [https://doi.org/10.1016/j.rse.2017.07.002](https://doi.org/10.1016/j.rse.2017.07.002). (ONLY for SCS+C function)
