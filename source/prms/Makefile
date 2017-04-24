# prms libraray makefile
# $Id: Makefile 8247 2016-04-14 22:53:07Z rsregan $

include ../makelist
PRMSLIB = $(LIBDIR)/libprmsiv.a
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

#
# Define all object files which make up the library
#
LIBOBJS = \
	basin.o \
        cascade.o \
        prms_time.o \
	climateflow.o \
        soltab.o \
        obs.o \
        climate_hru.o \
        temp_1sta_laps.o \
        temp_dist2.o \
        xyz_dist.o \
        ide_dist.o \
        precip_1sta_laps.o \
        precip_dist2.o \
        ddsolrad.o \
        ccsolrad.o \
        transp_tindex.o \
        transp_frost.o \
        frost_date.o \
        potet_jh.o \
        potet_hamon.o \
        potet_pan.o \
        potet_pm.o \
        potet_pm_sta.o \
        potet_pt.o \
        potet_hs.o \
        write_climate_hru.o \
        intcp.o \
        snowcomp.o \
        srunoff.o \
        soilzone.o \
        gwflow.o \
        routing.o \
        strmflow.o \
        strmflow_in_out.o \
        muskingum.o \
        utils_prms.o \
        basin_sum.o \
        nhru_summary.o \
        water_balance.o \
        map_results.o \
        subbasin.o

all: lib

install: lib

lib: $(PRMSLIB)

$(PRMSLIB): $(LIBOBJS)
	$(RM) $(PRMSLIB)
	$(AR) $(PRMSLIB) $(LIBOBJS)
	$(RANLIB) $(PRMSLIB)

clean:
	$(RM) $(PRMSLIB)
	$(RM) *.o *.mod *~

gwflow.o: gwflow.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_intcp.mod prms_srunoff.mod prms_soilzone.mod prms_cascade.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) gwflow.f90

utils_prms.o: utils_prms.f90 prms_module.mod prms_basin.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) utils_prms.f90

prms_time.o: prms_time.f90 prms_module.mod prms_basin.mod
	$(FC) -c $(FFLAGS) prms_time.f90

intcp.o: intcp.f90 prms_module.mod prms_climatevars.mod prms_basin.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) intcp.f90

map_results.o: map_results.f90 prms_module.mod prms_basin.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) map_results.f90

soltab.o: soltab.f90 prms_module.mod prms_basin.mod
	$(FC) -c $(FFLAGS) soltab.f90

precip_1sta_laps.o: precip_1sta_laps.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) precip_1sta_laps.f90

transp_tindex.o: transp_tindex.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) transp_tindex.f90

temp_1sta_laps.o: temp_1sta_laps.f90 prms_climatevars.mod prms_basin.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) temp_1sta_laps.f90

strmflow.o: strmflow.f90 prms_module.mod prms_basin.mod prms_flowvars.mod prms_gwflow.mod prms_srunoff.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) strmflow.f90

strmflow_in_out.o: strmflow_in_out.f90 prms_module.mod prms_basin.mod prms_obs.mod prms_flowvars.mod prms_gwflow.mod prms_set_time.mod prms_srunoff.mod
	$(FC) -c $(FFLAGS) strmflow_in_out.f90

muskingum.o: muskingum.f90 prms_module.mod prms_basin.mod prms_obs.mod prms_flowvars.mod prms_gwflow.mod prms_set_time.mod prms_routing.mod prms_srunoff.mod
	$(FC) -c $(FFLAGS) muskingum.f90

potet_jh.o: potet_jh.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) potet_jh.f90

potet_pm.o: potet_pm.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_climate_hru.mod prms_soltab.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) potet_pm.f90

potet_pm_sta.o: potet_pm_sta.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_climate_hru.mod prms_soltab.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) potet_pm_sta.f90

potet_pt.o: potet_pt.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) potet_pt.f90

potet_hs.o: potet_hs.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) potet_hs.f90

transp_frost.o: transp_frost.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) transp_frost.f90
 
frost_date.o: frost_date.f90 prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) frost_date.f90
 
potet_pan.o: potet_pan.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) potet_pan.f90

potet_hamon.o: potet_hamon.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) potet_hamon.f90

write_climate_hru.o: write_climate_hru.f90  prms_module.mod prms_climatevars.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) write_climate_hru.f90

water_balance.o: water_balance.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_soilzone.mod prms_srunoff.mod prms_gwflow.mod prms_set_time.mod prms_cascade.mod prms_intcp.mod prms_snow.mod
	$(FC) -c $(FFLAGS) water_balance.f90

climate_hru.o: climate_hru.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_soltab.mod
	$(FC) -c $(FFLAGS) climate_hru.f90

climateflow.o: climateflow.f90  prms_module.mod prms_basin.mod
	$(FC) -c $(FFLAGS) climateflow.f90

basin.o: basin.f90  prms_module.mod
	$(FC) -c $(FFLAGS) basin.f90

cascade.o: cascade.f90  prms_module.mod prms_basin.mod
	$(FC) -c $(FFLAGS) cascade.f90

obs.o: obs.f90  prms_module.mod prms_basin.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) obs.f90

temp_dist2.o: temp_dist2.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) temp_dist2.f90

precip_dist2.o: precip_dist2.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) precip_dist2.f90

ddsolrad.o: ddsolrad.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) ddsolrad.f90

ccsolrad.o: ccsolrad.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_soltab.mod prms_set_time.mod prms_obs.mod
	$(FC) -c $(FFLAGS) ccsolrad.f90

snowcomp.o: snowcomp.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_intcp.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) snowcomp.f90

srunoff.o: srunoff.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_intcp.mod prms_snow.mod prms_cascade.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) srunoff.f90

soilzone.o: soilzone.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_intcp.mod prms_cascade.mod prms_set_time.mod prms_snow.mod prms_srunoff.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) soilzone.f90

basin_sum.o: basin_sum.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_intcp.mod prms_snow.mod prms_srunoff.mod prms_gwflow.mod prms_set_time.mod prms_muskingum.mod
	$(FC) -c $(FFLAGS) basin_sum.f90

nhru_summary.o: nhru_summary.f90  prms_module.mod prms_basin.mod prms_set_time.mod
	$(FC) -c $(FFLAGS) nhru_summary.f90

routing.o: routing.f90  prms_module.mod prms_basin.mod prms_set_time.mod prms_climatevars.mod prms_set_time.mod prms_gwflow.mod
	$(FC) -c $(FFLAGS) routing.f90

subbasin.o: subbasin.f90  prms_module.mod prms_basin.mod prms_climatevars.mod prms_intcp.mod prms_snow.mod prms_srunoff.mod prms_set_time.mod prms_soilzone.mod prms_gwflow.mod
	$(FC) -c $(FFLAGS) subbasin.f90

ide_dist.o: ide_dist.f prms_module.mod prms_basin.mod prms_set_time.mod prms_climatevars.mod prms_obs.mod
	$(FC) -c $(FFLAGS) ide_dist.f

xyz_dist.o: xyz_dist.f prms_module.mod prms_basin.mod prms_set_time.mod prms_climatevars.mod prms_obs.mod
	$(FC) -c $(FFLAGS) xyz_dist.f

prms_climatevars.mod: climateflow.o
prms_flowvars.mod: climateflow.o
prms_gwflow.mod: gwflow.o
prms_snow.mod: snowcomp.o
prms_intcp.mod: intcp.o
prms_module.mod:
	$(CD) ../gsflow;make prms_module.mod
	$(CP) ../gsflow/prms_module.mod .
prms_obs.mod: obs.o
prms_basin.mod: basin.o
prms_soltab.mod: soltab.o
prms_muskingum.mod: muskingum.o
prms_set_time.mod: prms_time.o

