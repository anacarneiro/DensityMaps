# DensityMaps

### Mapping the global distribution of seabird populations: a framework for integrating tracking, demographic and phenological datasets.

This is a step-by-step guide to do the analyses presented in the manuscript. An example of tracking data is available for the Tristan Albatross. This example dataset can only be used to run the scripts below and cannot be used for dissemination activities (e.g. publication) or any other purpose without the approval of data owners. All datasets are available for request in the Seabird Tracking Database (http://seabirdtracking.org/). The scripts were developed assuming multiple species (i.e you can have multiple species files to run your analysis). Metadata on phenology and demography models are available for all species used in the manuscript. Example outputs from each step of the methods are available in the scripts_results folder (warning: the raster files are very large).

## Steps:
1. **Download data from the Seabird Tracking Database**
1. **Create one csv file per species**
1. **Follow the instructions and run the R files** (in the corresponding order):
    * 01_demography 
    * 02_cleaning_data
    * 03_kernels
1. **Do a bootstrap analysis** – following instructions as in Lascelles et al. (2016) and Oppel et al. (2018)
1. **Follow the instructions and run the R files** (in the corresponding order):
    * 04_ combining_selecting_renaming
    * 05_land_mask
1. **Metadata creation based on phenology**
    * Fill WORKING_monthly_datatypes (inside folder metadata_files -> 06_metadata_phenology) for your species based on phenology. This is only for illustrative purposes to help filling the next step. The idea of the template is similar for most species, with differences in the format for biennial and annual species. 
    * Fill metadata_monthly_datatypes (inside folder metadata files -> metadata_phenology) for your species. Use the template in sheet 1 to fill out the metadata structure using phenology (average start dates and number of days) from the WORKING_monthly_datatypes table. Read into R, and follow the next steps.
    * **Alternatively**, use the Shiny app available at https://github.com/lizziepear/seabird-phenology to explore the phenology metadata tables and download the full phenology metadata table for your population as a csv (this enables you to skip script 06, except for adding in the 'value' field from the results of your demography model).
1. **Follow the instructions and run the R files** (in the corresponding order):
    * 06_combining_metadata_files
    * 07_monthly_equations
    * 08_quarter_combination
    * 09_sum_demClasses
    * 10_aggregate_5by5_grid
    * 11_year_combination

## License Information:
* The **R scripts** are licensed under the GNU General Public License version 3 (see LICENSE.md).
* The **example tracking data for Tristan Albatross** are licensed under **'all rights reserved'**. You may **only** use this data for personal use to run the scripts provided. You may **not** use the data for dissemination activities (for example, in published papers or showing/storing the data in an online platform). If you wish to use these data for any other purposes, please contact the authors of the paper, or make a direct request for the data via the Seabird Tracking Database.
