# Example-Raman-spectral-analysis

Can be used as example for a data analysis workflow as suggested by Guo et al "Tutorial on Chemometric Analysis in Raman Spectroscopy: from Experimental Design to Machine Learning based Modelling". It also generates all the images related to the mice dataset in the publication.

While the Script "Analysis_code.R" contains the analysis code for thze example data, the file "function_definitions.R" contains the functions definitions.

The files "DATA_dp_wc_bc.csv" and "DATA_dp_wc_bc.RData" contain Raman spectra from mice tissues, averaged from each Raman map after despiking, wavenumber calibration and baseline correction.  Details of these pre-processing steps are referred to S. Guo, et al., Journal of Chemometrics, 2020, 34: e3202. Intensities at different wavenumbers are saved column-wise. The ID, annotation, and type of the tissues are included as separate columns in the file.
