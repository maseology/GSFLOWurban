# GSFLOW Urban

An urban development water balance analysis expansion to [GSFLOW: coupled groundwater and surface-water flow model](https://water.usgs.gov/ogw/gsflow/).  

Built upon (and backward compatible to) GSFLOW version 1.2.2.

**Model capabilities:**

 * Division of HRU into three independent water balance features:  
   1. _Connected_ impervious areas representative of on-grade impervious areas
   2. _Disconnected_ impervious areas removed from the runon process meant to represent elevated structures such as roof tops
   3. Pervious areas
 * Conceptual _"Sewershed"_ independent of overland and groundwater pathways used to simulate urban storm water drainage
 * An infiltration storage feature used to investigate sustainable development strategies
 * A variety of logical pathways needed to routed water among various storage reservoirs
 * Groundwater interaction with:  
   1. Infiltration storage reservoirs in order to quantify long-term storage/retention potential
   2. Sewershed elements to simulation groundwater infiltration (i.e., storm sewer I&I)

![GSFLOW urban flow pathways](/doc/pathways_170515.png)

### Input instructions
...can be found [here](/doc/input_instructions.pdf).

### Current (Beta) version 0.3
**Task list:**

 - [x] Write & compile code
 - [x] Preliminary code testing (check for water balance closure)
 - [x] Create model input instructions
 - [ ] Complete model manual
 - [ ] Rigorous code (crash) testing
 - [ ] Sample problem
 
### Code modifications

The following table lists the model files that have been altered/modified as part of the GSFLOW urban extension. Code changes to existing USGS model files have been commented by the contributors' initials.  

Not all changes were implemented for the creation of GSFLOW urban, but modifications to the code were added to improve the legibility of simulation warning messages, such as increasing the digits when writing to the console

Users are free to recompile their own version of GSFLOW urban by simply replacing their GSFLOW project files with the files listed below. The contributors have assured that all changes made to the GSFLOW code will not affect the original code implementation, meaning the GSFLOW urban will remain backward compatible such that any previously-built GSFLOW model will continue to run using GSFLOW urban.

Original model file | Modified model file
------------------- | -------------------
gsflow/gsflow_budget.f90 | gsflow_budget_gu.f90
gsflow/gsflow_modflow.f | gsflow_modflow_gu.f
gsflow/gsflow_prms.f90 | gsflow_prms_gu.f90
gsflow/gsflow_prms2mf.f90 | gsflow_prms2mf_gu.f90
modflow/de47_NWT.f | de47_NWT_gu.f
modflow/gwf2bas7_NWT.f | gwf2bas7_NWT_gu.f
modflow/gwf2uzf1_NWT.f | gwf2uzf1_NWT_gu.f
prms/gwflow.f90 | gwflow_gu.f90
prms/map_results.f90 | map_results_gu.f90
prms/soilzone.f90 | soilzone_gu.f90
prms/srunoff.f90 | srunoff_gu.f90
prms/water_balance.f90 | water_balance_gu.f90
_none_ | srunoff_urban_gu.f90


### License

GSFLOWurban hosted on GitHub is released under the MIT license requiring the preservation of copyright and license notices. Licensed works, modifications, and larger works may be distributed under different terms and without source code. This code extension follows strict adherence to the USGS [Software User Rights Notice](/USGS_Software_User_Rights_Notice.md), which can also be found [here](https://water.usgs.gov/software/help/notice/).

### Contributors

Mason Marchidon P.Eng M.A.Sc, Hydrologist for the [Oak Ridges Moraine Groundwater Program](http://oakridgeswater.ca/)

Copyright © 2017

### References and Acknowledgments

Markstrom, S.L., Niswonger, R.G., Regan, R.S., Prudic, D.E., and Barlow, P.M., 2008, [GSFLOW-Coupled Ground-water and Surface-water FLOW model based on the integration of the Precipitation-Runoff Modeling System (PRMS) and the Modular Ground-Water Flow Model (MODFLOW-2005): U.S. Geological Survey Techniques and Methods 6-D1, 240 p.](https://pubs.usgs.gov/tm/tm6d1/)

A special thanks to the [Lake Simcoe Region Conservation Authority](http://www.lsrca.on.ca/) for their start-up funding contribution.
