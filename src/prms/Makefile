# PRMS library makefile

include ../makelist
PRMSLIB = $(LIBDIR)/libprms.a
MMFLIB	 = $(LIBDIR)/libmmf.a
LIBS	= $(PRMSLIB) $(MMFLIB) $(FLIBS)

####################################################
# Rules for targets
####################################################
all: $(PRMSLIB)

$(TARGET): $(PRMSLIB)
	$(RM) $(TARGET)

$(MMFLIB):
	$(CD) $(MMFDIR);make

####################################################
# Rules for targets
####################################################

#
# Define all object files which make up the library
#


LIBOBJS = \
        prms_constants_gu.o \
        basin.o \
        climateflow_gu.o \
        cascade_gu.o \
        soltab.o \
        convert_params.o \
        prms_time_gu.o \
        obs_gu.o \
        climate_hru.o \
        potet_jh.o \
        potet_pt.o \
        potet_hs.o \
        potet_pm.o \
        potet_pm_sta.o \
        potet_pan.o \
        potet_hamon.o \
        ddsolrad.o \
        ccsolrad.o \
        ide_dist.o \
        xyz_dist.o \
        precip_1sta_laps.o \
        temp_map.o \
        precip_map.o \
        precip_dist2.o \
        temp_1sta_laps.o \
        temp_dist2.o \
        transp_frost.o \
        transp_tindex.o \
        frost_date.o \
        glacr_melt.o \
        intcp.o \
        snowcomp_gu.o \
        srunoff_gu.o \
        srunoff_urban_gu.o \
        soilzone_gu.o \
        gwflow_gu.o \
        water_use_read.o \
        dynamic_param_read.o \
        water_balance_gu.o \
        routing.o \
        strmflow.o \
        strmflow_in_out.o \
        muskingum.o \
        muskingum_lake.o \
        subbasin.o \
        stream_temp.o \
        map_results_gu.o \
        nhru_summary.o \
        nsub_summary.o \
        nsegment_summary.o \
        basin_summary.o \
        write_climate_hru.o \
        prms_summary.o \
        basin_sum.o \
        utils_prms.o

install: lib

lib: $(PRMSLIB)

$(PRMSLIB): $(LIBOBJS)
	$(RM) $(PRMSLIB)
	$(AR) $(PRMSLIB) $(LIBOBJS)
	$(RANLIB) $(PRMSLIB)

clean:
	$(RM) *.o *.mod *~

#
# Define all object files which make up the library
#

prms_constants_gu.o: prms_constants_gu.f90
	$(FC) -c $(FFLAGS) prms_constants_gu.f90

basin_sum.o: basin_sum.f90 prms_module.mod prms_flowvars.mod prms_intcp.mod prms_snow.mod prms_srunoff.mod prms_gwflow.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_muskingum.mod prms_constants.mod
	$(FC) -c $(FFLAGS) basin_sum.f90

subbasin.o: subbasin.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_set_time.mod prms_intcp.mod prms_srunoff.mod prms_soilzone.mod prms_gwflow.mod prms_snow.mod prms_climatevars.mod prms_muskingum_lake.mod prms_constants.mod prms_climate_hru.mod
	$(FC) -c $(FFLAGS) subbasin.f90

ddsolrad.o: ddsolrad.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) ddsolrad.f90

ccsolrad.o: ccsolrad.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) ccsolrad.f90

snowcomp_gu.o: snowcomp_gu.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_climatevars.mod prms_set_time.mod prms_intcp.mod prms_constants.mod prms_soltab.mod prms_climate_hru.mod
	$(FC) -c $(FFLAGS) snowcomp_gu.f90

gwflow_gu.o: gwflow_gu.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_intcp.mod prms_srunoff.mod prms_soilzone.mod prms_cascade.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) gwflow_gu.f90

utils_prms.o: utils_prms.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) utils_prms.f90

prms_summary.o: prms_summary.f90 prms_module.mod prms_climatevars.mod prms_flowvars.mod prms_set_time.mod prms_obs.mod prms_intcp.mod prms_snow.mod prms_srunoff.mod prms_soilzone.mod prms_gwflow.mod prms_constants.mod
	$(FC) -c $(FFLAGS) prms_summary.f90

muskingum.o: muskingum.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_set_time.mod prms_obs.mod prms_srunoff.mod prms_gwflow.mod prms_routing.mod  prms_constants.mod
	$(FC) -c $(FFLAGS) muskingum.f90

intcp.o: intcp.f90 prms_module.mod prms_basin.mod prms_obs.mod prms_climatevars.mod prms_flowvars.mod prms_set_time.mod prms_water_use.mod prms_constants.mod
	$(FC) -c $(FFLAGS) intcp.f90

map_results_gu.o: map_results_gu.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) map_results_gu.f90

nhru_summary.o: nhru_summary.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) nhru_summary.f90

nsub_summary.o: nsub_summary.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) nsub_summary.f90

nsegment_summary.o: nsegment_summary.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) nsegment_summary.f90

basin_summary.o: basin_summary.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) basin_summary.f90

soltab.o: soltab.f90 prms_module.mod prms_basin.mod prms_constants.mod
	$(FC) -c $(FFLAGS) soltab.f90

frost_date.o: frost_date.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) frost_date.f90

precip_1sta_laps.o: precip_1sta_laps.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) precip_1sta_laps.f90

temp_map.o: temp_map.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) temp_map.f90

precip_map.o: precip_map.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) precip_map.f90

transp_tindex.o: transp_tindex.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) transp_tindex.f90

transp_frost.o: transp_frost.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) transp_frost.f90

temp_1sta_laps.o: temp_1sta_laps.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) temp_1sta_laps.f90

temp_dist2.o: temp_dist2.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) temp_dist2.f90

precip_dist2.o: precip_dist2.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) precip_dist2.f90

strmflow.o: strmflow.f90 prms_module.mod prms_basin.mod prms_obs.mod prms_flowvars.mod prms_gwflow.mod prms_srunoff.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) strmflow.f90

strmflow_in_out.o: strmflow_in_out.f90 prms_module.mod prms_basin.mod prms_obs.mod prms_flowvars.mod prms_gwflow.mod prms_srunoff.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) strmflow_in_out.f90

stream_temp.o: stream_temp.f90 prms_module.mod prms_basin.mod prms_obs.mod prms_flowvars.mod prms_routing.mod prms_set_time.mod prms_snow.mod prms_soltab.mod prms_constants.mod
	$(FC) -c $(FFLAGS) stream_temp.f90

potet_jh.o: potet_jh.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_jh.f90

potet_pt.o: potet_pt.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_climate_hru.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_pt.f90

potet_hs.o: potet_hs.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_hs.f90

potet_pm.o: potet_pm.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_climate_hru.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_pm.f90

potet_pm_sta.o: potet_pm_sta.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_climate_hru.mod prms_soltab.mod prms_set_time.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_pm_sta.f90

potet_pan.o: potet_pan.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_obs.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_pan.f90

potet_hamon.o: potet_hamon.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) potet_hamon.f90

write_climate_hru.o: write_climate_hru.f90 prms_module.mod prms_set_time.mod prms_climatevars.mod prms_constants.mod
	$(FC) -c $(FFLAGS) write_climate_hru.f90

climate_hru.o: climate_hru.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) climate_hru.f90

cascade_gu.o: cascade_gu.f90 prms_module.mod prms_basin.mod prms_constants.mod
	$(FC) -c $(FFLAGS) cascade_gu.f90

basin.o: basin.f90 prms_module.mod prms_constants.mod
	$(FC) -c $(FFLAGS) basin.f90

obs_gu.o: obs_gu.f90 prms_module.mod prms_set_time.mod prms_basin.mod prms_constants.mod
	$(FC) -c $(FFLAGS) obs_gu.f90

srunoff_gu.o: srunoff_gu.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_climatevars.mod prms_cascade.mod prms_intcp.mod prms_snow.mod prms_set_time.mod prms_constants.mod prms_climate_hru.mod srunoff_urban.mod
	$(FC) -c $(FFLAGS) srunoff_gu.f90

srunoff_urban_gu.o: srunoff_urban_gu.f90 prms_module.mod prms_basin.mod prms_constants.mod
	$(FC) -c $(FFLAGS) srunoff_urban_gu.f90	

climateflow_gu.o: climateflow_gu.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_constants.mod
	$(FC) -c $(FFLAGS) climateflow_gu.f90

routing.o: routing.f90 prms_module.mod prms_basin.mod prms_gwflow.mod prms_flowvars.mod prms_set_time.mod prms_srunoff.mod prms_climatevars.mod prms_water_use.mod prms_constants.mod
	$(FC) -c $(FFLAGS) routing.f90

prms_time_gu.o: prms_time_gu.f90 prms_module.mod prms_basin.mod prms_constants.mod
	$(FC) -c $(FFLAGS) prms_time_gu.f90

water_balance_gu.o: water_balance_gu.f90 prms_module.mod prms_basin.mod prms_srunoff.mod prms_flowvars.mod prms_gwflow.mod prms_climatevars.mod prms_set_time.mod prms_cascade.mod prms_intcp.mod prms_snow.mod prms_soilzone.mod prms_constants.mod srunoff_urban.mod
	$(FC) -c $(FFLAGS) water_balance_gu.f90

ide_dist.o: ide_dist.f prms_module.mod prms_basin.mod prms_set_time.mod prms_climatevars.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) ide_dist.f

xyz_dist.o: xyz_dist.f prms_module.mod prms_basin.mod prms_set_time.mod prms_climatevars.mod prms_obs.mod prms_constants.mod
	$(FC) -c $(FFLAGS) xyz_dist.f

muskingum_lake.o: muskingum_lake.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_set_time.mod prms_obs.mod prms_routing.mod prms_srunoff.mod prms_gwflow.mod prms_soilzone.mod prms_constants.mod
	$(FC) -c $(FFLAGS) muskingum_lake.f90

dynamic_param_read.o: dynamic_param_read.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_climatevars.mod prms_flowvars.mod prms_potet_jh.mod prms_potet_pm.mod prms_potet_hs.mod prms_potet_pt.mod prms_potet_hamon.mod transp_tindex.o transp_frost.o prms_intcp.mod prms_snow.mod prms_srunoff.mod prms_soilzone.mod prms_climate_hru.mod prms_constants.mod
	$(FC) -c $(FFLAGS) dynamic_param_read.f90

water_use_read.o: water_use_read.f90 prms_module.mod prms_basin.mod prms_set_time.mod prms_flowvars.mod prms_constants.mod
	$(FC) -c $(FFLAGS) water_use_read.f90

convert_params.o: convert_params.f90 prms_module.mod
	$(FC) -c $(FFLAGS) convert_params.f90

glacr_melt.o: glacr_melt.f90 prms_snow.mod prms_intcp.mod prms_soltab.mod prms_constants.mod
	$(FC) -c $(FFLAGS) glacr_melt.f90

prms_climatevars.mod: climateflow_gu.o
prms_flowvars.mod: climateflow_gu.o
prms_gwflow.mod: gwflow_gu.o
prms_module.mod:
	$(CD) ../gsflow;make prms_module.mod
	$(CP) ../gsflow/prms_module.mod .
prms_obs.mod: obs_gu.o
prms_basin.mod: basin.o
prms_soltab.mod: soltab.o
prms_muskingum.mod: muskingum.o
prms_intcp.mod: intcp.o
prms_snow.mod: snowcomp_gu.o
prms_cascade.mod: cascade_gu.o
prms_srunoff.mod: srunoff_gu.o
srunoff_urban.mod: srunoff_urban_gu.o
prms_soilzone.mod: soilzone_gu.o
prms_routing.mod: routing.o
prms_water_use.mod: water_use_read.o
prms_set_time.mod: prms_time_gu.o
prms_glacr.mod: glacr_melt.o
prms_constants.mod: prms_constants_gu.o
prms_climate_hru.mod: climate_hru.o
prms_potet_jh.mod: potet_jh.o
prms_potet_pm.mod: potet_pm.o
prms_potet_pt.mod: potet_pt.o
prms_potet_hs.mod: potet_hs.o
prms_potet_hamon.mod: potet_hamon.o
gsfmodflow.mod: $(GSFLOWDIR)/gsfmodflow.mod
	$(RM) gsfmodflow.mod
	$(LN) $(GSFLOWDIR)/gsfmodflow.mod gsfmodflow.mod
gwfagmodule.mod: $(MODFLOWDIR)/gwfagmodule.mod
	$(RM) gwfagmodule.mod
	$(LN) $(MODFLOWDIR)/gwfagmodule.mod gwfagmodule.mod
	
