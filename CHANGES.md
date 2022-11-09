# Changelog
All notable changes to this project will be documented in this file

## [Unreleased]

- New build system based on autoconf
- Major bug fix for relative azimuth angle and scattering angle
- New Sensors
    - JPSS-2 / NOAA-21
    - VGAC Fusion
    

## [v1.0.0]

- UMD surface emission option fixed
- Shadow mask disabled
- Greatly increased processing speed
- Added support for GOES-18
- Added RTTOV support for GOES-8
- Set up continuous integration (CI)
- ECM2 bug fixes: stddev precision and bin fix
- Catch SIGINT / SIGTERM and clean-up tmpdirs before exiting
- Better repeatibility (initialize memory)
- Added ISCCP-NG L1g support
- GEO subsets are cropped to latitude and longitude
- Improved time (changes solar geometry and reflectance)
- Updated cloud type algorithm for "Mixed" type
- Changes to ACHA quality flag
- Clip ACHA beta between 1 and 2
- Uncertainty in tau added to ACHA estimation
- Updated cloud phase / type post-ECM2 filtering
- Misc bug fixes
- Added support for FY3D


