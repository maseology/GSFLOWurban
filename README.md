# GSFLOWurban

An urban development water balance analysis expansion to [GSFLOW: coupled groundwater and surface-water flow model](https://www.usgs.gov/software/gsflow-coupled-groundwater-and-surface-water-flow-model). Built upon (and backward compatible to) [GSFLOW version 2.2.1](https://water.usgs.gov/water-resources/software/gsflow/GSFLOW_Release_Notes_2.2.1.pdf). 

[See *Draft* Manual addendum](/doc/GSFU_man_Jan19.pdf).

## Model capabilities:

 * Division of HRU into three independent water balance features:  
   1. _Connected_ impervious areas representative of on-grade impervious areas
   2. _Disconnected_ impervious areas removed from the runon process meant to represent elevated structures such as roof tops
   3. Pervious areas
 * Conceptual _"Sewershed"_ independent of overland and groundwater pathways used to simulate urban storm water drainage
 * An infiltration storage feature used to investigate sustainable development strategies
 * A variety of logical pathways needed to route water among various storage reservoirs
 * Groundwater interaction with:  
   1. Infiltration storage reservoirs in order to quantify long-term storage/retention potential
   2. Sewershed elements to simulation groundwater infiltration (i.e., storm sewer I&I)

![GSFLOW urban flow pathways](/doc/pathways_181211.png)

## Input instructions
...can be found [here](/doc/input_instructions.pdf).

## Current (Beta) version 0.6.1

**Task list:**

 - [x] Write & compile code
 - [x] Preliminary code testing (check for water balance closure)
 - [x] Create model input instructions
 - [ ] Complete model manual
 - [ ] Rigorous code (crash) testing
 - [ ] Sample problem

## Additional processes added to GSFLOW

* SCS-CN (1972) runoff generation mechanism as applied to the SWAT model (Neitsch et.al., 2011) 
* Variable (sub-daily) time step Green and Ampt (1911) solution of Chu (1978)
* TOPMODEL (Beven et.al., 1995) as an alternative to cascading groundwater module in PRMS-only mode
 
## Code modifications

The following table lists the model files that have been altered/modified as part of the GSFLOWurban extension. Code changes to existing USGS model files have been commented by the contributors' initials.  

Not all changes were implemented for the creation of GSFLOW urban, but modifications to the code were added to improve the legibility of simulation warning messages, such as increasing the digits when writing to the console

Users are free to recompile their own version of GSFLOWurban by simply replacing their GSFLOW project files with the files listed below. The contributors have assured that all changes made to the GSFLOW code will not affect the original code implementation, meaning the GSFLOWurban will remain backward compatible such that any previously-built GSFLOW model will continue to run using GSFLOWurban.

Original model file | Modified model file
------------------- | -------------------
gsflow/gsflow_budget.f90 | gsflow_budget_gu.f90
gsflow/gsflow_modflow.f | gsflow_modflow_gu.f
gsflow/gsflow_prms.f90 | gsflow_prms_gu.f90
gsflow/gsflow_prms2mf.f90 | gsflow_prms2mf_gu.f90
gsflow/gsflow_sum.f90 | gsflow_sum_gu.f90
mmf/read_line.c | read_line_gu.c
mmf/setup_cont.c | setup_cont_gu.c
modflow/de47_NWT.f | de47_NWT_gu.f
modflow/gwf2bas7_NWT.f | gwf2bas7_NWT_gu.f
modflow/gwf2uzf1_NWT.f | gwf2uzf1_NWT_gu.f
modflow/NWT1_solver.f | NWT1_solver_gu.f
prms/cascade.f90 | cascade_gu.f90
prms/climateflow.f90 | climateflow_gu.f90
prms/gwflow.f90 | gwflow_gu.f90
prms/map_results.f90 | map_results_gu.f90
prms/obs.f90 | obs_gu.f90
prms/prms_constants.f90 | prms_constants_gu.f90
prms/prms_time.f90 | prms_time_gu.f90
prms/snowcomp.f90 | snowcomp_gu.f90
prms/soilzone.f90 | soilzone_gu.f90
prms/srunoff.f90 | srunoff_gu.f90
prms/water_balance.f90 | water_balance_gu.f90
_none_ | srunoff_urban_gu.f90

### Changes to parameter input ranges

A number of parameter ranges have been changed to satisfy specific needs. Setting parameter values beyond the default ranges set in Markstrom et.al. (2008) should be made with caution.

* *Imperv_stor_max* range increased to [0.0 to 5.0] inches
* *Soil_rechr_max* range increased to [0.00001 to 10.0] inches
* *Ssr2gw_rate* range increased from [0.0 to 1.0] to [0.0 to 10,000.0]; however, *Ssr2gw_rate* > 1.0 will only be accepted if *Ssr2gw_exp* is set to zero (=0.0). This modification allows the potential gravity drainage from the soil zone to be set explicitly as a rate and not as a function of gravity reservoir storage (see equation 59 in Markstrom et.al., 2008).
* *Sat_threshold* lower limit of range changed to [0.0 to 999.0] inches
* *Cecn_coef* range increased to [0.01 to 60.0] calories per degree Celsius above 0
* *ncol* range increased to [1 to 500,000]

# License

GSFLOWurban hosted on GitHub is released under the MIT license requiring the preservation of copyright and license notices. Licensed works, modifications, and larger works may be distributed under different terms and without source code. This code extension follows strict adherence to the USGS [Software User Rights Notice](/USGS_Software_User_Rights_Notice.md), which can also be found [here](https://water.usgs.gov/software/help/notice/).

# Contributors

Mason Marchildon P.Eng M.A.Sc, Hydrologist for the [Oak Ridges Moraine Groundwater Program](http://oakridgeswater.ca/)

Peter John Thompson P.Eng M.A.Sc, Senior Hydrologist, [GeoProcess Research Associates](https://geoprocess.com/)

Copyright © 2017

# References and Acknowledgements

Beven, K.J., R. Lamb, P.F. Quinn, R. Romanowicz, and J. Freer, 1995. TOPMODEL. In Singh V.P. editor, Computer Models of Watershed Hydrology. Water Resources Publications, Highland Ranch, CO: pp. 627–668.

Chu, S.T., 1978. Infiltration During an Unsteady Rain. Water Research Research 14(3). pp.461-466.

Green, W.H., G.A. Ampt, 1911. Studies of Soil Physics, 1: The Flow of Air and Water Through Soils. Journal of Agricultural Science 4(1). pp.1-24.

Markstrom, S.L., Niswonger, R.G., Regan, R.S., Prudic, D.E., and Barlow, P.M., 2008, [GSFLOW-Coupled Ground-water and Surface-water FLOW model based on the integration of the Precipitation-Runoff Modeling System (PRMS) and the Modular Ground-Water Flow Model (MODFLOW-2005): U.S. Geological Survey Techniques and Methods 6-D1, 240 p.](https://pubs.usgs.gov/tm/tm6d1/)

Neitsch, S.L., J.G. Arnold, J.R. Kiniry, J.R. Williams, 2011. Soil and Water Assessment Tool: Theoretical Documentation, version 2009: TR-406. 647pp.

United States Department of Agriculture: Soil Conservation Service, 1972. National Engineering Handbook, Section 4–Hydrology. USDA-SCS, Washington, D.C.

<br>

A special thanks to the [Lake Simcoe Region Conservation Authority](http://www.lsrca.on.ca/) for their start-up funding contribution.


# Release info

#### v0.6.1 February, 2019
* re-organized and tested sub-daily Green and Ampt code

#### v0.6 January, 2019
* added sub-daily precipitation input (thanks to @PJ-Thompson). *update to manual coming shortly*

#### v0.5 January, 2019
* changes to overall srunoff module: now urban hydrology mode can be combined with any srunoff module (i.e., srunoff_carea, srunoff_smidx, srunoff_scscn, srunoff_grnampt, etc.) 
* changes to SCS-CN module to reflect manual methodology
* include the PRMS-TOPMODEL option (this module had somehow been forgotten, oops **;)**

#### v0.4 December, 2018

* released draft manual
* changes to Green and Ampt module to allows the soil moisture accounting dictate infiltrability
* configured makefiles for compiling in Linux
* bug fixes and input range adjustments
* (began including release info)