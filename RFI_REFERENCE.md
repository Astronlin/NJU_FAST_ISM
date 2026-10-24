# Published FAST RFI references

The package catalog combines empirical observation profiles, published historical statistics and published harmonic patterns. Published frequency statistics are not automatic masks. Existing IC443 seeds and fitting settings are unchanged.

`data/fast_rfi_statistics_zhang2022.json` preserves all 5001 numeric rows from ScienceDB V2, CC0. Original `rfi_statis.txt`: 275335 bytes; MD5 514051d82b50c92bdef040ba206c5826. Sources and retrieval checksums are recorded in JSON. The original DOI redirects to 10.57760/sciencedb.j00001.00313.

These statistics cover July 29–September 11, 2021, approximately 300 hours, with 0.1 MHz sampling and an approximately 0.324 K threshold. They include Galactic H I near 1420.5 MHz, explicitly identified in the paper; never infer RFI from high occupancy near the astronomical line. Stored intensity zeros are not proof of RFI absence. The 2021 January/February IC443 observation predates the statistics.

The catalog records nearest-bin comparisons for every IC443 seed. A nonzero tabulated intensity in the same frequency bin is supporting historical co-location, not proof of the same emitter. Published external-service bands from the observatory's July 2021 tower-site report do not overlap the IC443 band. The report used a separate monitoring antenna and cannot inventory all receiver-internal interference.

HiFAST III section 6.1 reports the 8.1/16.2 MHz compressor pattern before late July 2021, when shielding removed it. This is consistent with IC443 spacings and its date, but does not establish every peak's origin. The 0.5 MHz ADC pattern is beam/polarization dependent; the 0.37 MHz L-band occurrence was reported for a few days in late November 2022. Neither supplies universal absolute masking frequencies.

Use `load_published_statistics(catalog, base_dir)` and `lookup_published_statistics(frequency_MHz, catalog, base_dir)` from `fastfish.rfi_catalog` for advisory lookups. Keep the observation's protected H I interval and require independent spectral/time/beam evidence before masking. No existing run products are rewritten by a reference update.

Detection policy: calibration now searches both observation-profile seeds and historical bins with reported positive flagged intensities. A mask requires an observed positive highpass peak above max(0.15 K-equivalent, 6 robust scatter estimates) in the short time-chunk spectrum. Published seeds search within 0.1 MHz; empirical seeds allow 0.8 MHz drift. Absent/weak priors are recorded as nondetections and never create a fallback mask. The protected H I interval is excluded from the search, and padding is clipped at its boundary. Short-time/128-channel detection can miss narrow, transient, weak or uncatalogued RFI; masks remain candidate classifications, not guarantees. Existing products/masks are unchanged until processing is explicitly rerun.
