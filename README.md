# **CL**ouds from **AV**h**R**r (and more) e**X**tended (CLAVR-x)

## Directory Structure

  build/ : contains things which pertain to compiling and building the package

  docs/ : contains various relevant documentation

  src/ : contains the package source code files
    acha/ : awg cloud height algorithm codes  
    asos/ : automatic surface observing system codes  
    baseline_cloud_height/ : baseline cloud height and 11um emissivity  
    baseline_cloud_mask/ : GOES-R baseline cloud mask  
    ccl/ : cloud cover layers  
    cloud_base/ : cloud base algorithm codes  
    cloud_mask/ : enterprise cloud mask codes  
    cloud_type/ : clavrx cloud typing codes  
    cx_dncomp/ : the day / night cloud optical and microphysical codes  
    cx_sds_io/ : general I/O code for HDF4-, HDF5-, and NetCDF- format datasets  
    dark_composite/ : code for geostationary dark-sky vis composites  
    main/ : houses most of the CLAVR-x processing system code  
    misc/ : miscellaneous code that does not seem to fit elsewhere  
    muri/ : aerosol retrieval code  
    pfaast/ : fast IR RTM (Pressure-Layer Fast Algo. for Atmos. Transmittance)  
    rttov/ : code to interface with the RTTOV RTM (Radiative Transfer for TOVS)  
    sasrab/ : the solar insolation codes  

  tests/ : contains a variety of tests that can be run manually or during CI/CD  
  Requres py.test; `test.sh` runs `py.test` in the `tests/` directory.  
  Sample output from select tests are saved for each release at  
  `/ships19/cloud/archive/clavrx_test_data/version_granules`  

  run/ : will be present after a successful build, and contains the built  
  executables and example/template run configuration files  
