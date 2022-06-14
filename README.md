# **CL**ouds from **AV**h**R**r (and more) e**X**tended (CLAVR-x)

## Directory Structure

    acha: awg cloud height algorithm codes
    asos: automatic surface observing system codes
    baseline_cloud_mask: GOES-R baseline cloud mask 
    ccl: cloud cover layers
    clavrx_bin: directory that holds the clavrx executables
    cloud_base: cloud base algorithm codes
    cloud_mask: enterprise cloud mask codes
    cloud_type: clavrx cloud typing codes
    main_src: houses the clavr-x processing system code and this where the Makefile exists
    pfaast: the fast IR RTM codes
    sasrab: the solar insolation codes
    dncomp: the day / night cloud optical and microphysical codes.  Note dncomp is compiled outside of the other CLAVR-x codes
    dark_composite: code that makes a seperate executable for make geostationary dark-sky vis composites

## Other Files

    clavrx_options_example: an example template of the clavrx options file needed by clavrxorb
    file_list_example: an example file list file needed by clavrxorb


## Testing

Requres py.test

`test.sh` simply runs `py.test` in the `test/` directory.

Sample output from select tests are saved for each release at `/ships19/cloud/archive/clavrx_test_data/version_granules`

