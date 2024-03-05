
<b>WYTHAM GREAT TITS - animal models</b>
==============================

Analysis for animal model paper...

<b>To use this repository:</b>
- In the terminal navigate to the folder where you want to install the repository. Then type `git clone https://github.com/carysvjones/AnimalMod_Envir.git`
- Or just download as .zip file, using the drop down in `<> Code`

<br>

<b>Run main analysis:</b>

1. Open the `AnimalMod_Envir.Rproj` file.
   
2. First open and run `analysis/0_pre_prep_data.R` - will output clean data files that are needed to run the next script,
you only have to run this once to create the datasets you need from the raw data (I just included it here for transparency so you can see what I did to the data).
Uses functions from `R/clean_data.R`, for cleaning up the data, mostly removing or correcting individual breeding attempts etc

3. Then once you've got clean data open and run `analysis/1_Prep_data.Rmd` - it preps the data for running animal models. Also run `analysis/1.1_Create_matrices.R` to create spatial and breeding environment similarity matrices for models.

4. Next run code in `analysis/2_Run_models.Rmd` to run the animal models using Asreml-R, this requires a licence, and models take a little while to run. Output from the models is already in the `model_output` folder, so the next scripts looking at results can still be run.

5. `analysis/3_Model_output.Rmd` goes through extracting variance componenets from all models and merging to tables, as well as plotting the output.

<br>

<b>Supplementary Materials:</b>

6. `analysis/3.1_Regression_analysis.R` is analysis for Section 1 of the Supplementary Info, including data analysis and plots.
7. `analysis/3.2_Supplementary.R` creates all plots in Section 2 of the Supplementary Info.
