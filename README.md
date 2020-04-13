# Cytof analysis Report

The aim of this Rmarkdown code is to provide the starting point for your CyTOF analysis. 

For more documentation and to view an example final report please view the ![documentation](https://acribbs.github.io/cytof_report/)

The report uses the Bioconductor packaged Flowcore and CATALYST to prodivde summary statistics and basic high dimensional cluserting and visualisations. These can be useful for determining the quality of your CyTOF data. Once the basic report has been generated you are then free to extend the report by adding your own project specific analysis code. 

## Configuring the report

The report requires you to modify the values within the config.yml file to customise the running of the workflow and add the required information within the *design_*.csv* file that specifies the sample name and meta data associated with your project.

FCS files are then placed within the data/ directory and the report ran using the Rstudio and the "Build Website" functionality.

# Step by step

* Clone the repository:

`git clone https://github.com/Acribbs/cytof_report.git`

* Navigate to the directory and rename the cytof_report.Rproj file to something of your choosing and then double click it to open Rstudio.

* Modify the design_*.csv and config.yml files to reflect your CytOF samples. You can have multiple design_* files in the folder.

* Once Rstudio is open then configure the build tools to run a website and then click "Build Website"

![Location of Build Website in Rstudio](https://raw.githubusercontent.com/Acribbs/deseq2_report/master/img/build_img.png)

Make sure rmarkdown is installed in your library and hit the build tab on the environment window and then click "Build Website". When the website has finished building a window will pop up with the rendered site. The final report is in the directory "Final_report" and can be accessed by opening the index.html file in a web browser.
