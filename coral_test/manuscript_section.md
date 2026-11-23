# Long-record test of the bidecadal line with coral archives

*Draft section for the manuscript. It fits a version of the paper that states the solar link as a hypothesis.*

### B.1 Motivation

The instrumental record (ERSSTv6, 1850–2025, N = 176 yr) shows a weak spectral line near 20 yr (amplitude about
0.03 K). Two periods are resolvable only if the record length N exceeds 1/|1/T1 − 1/T2| (Rayleigh criterion; with
noise, somewhat longer records are needed): separating 18.6 yr (lunar nodal cycle) from 20 yr requires about
270 yr, 20 from 22 yr (Hale cycle) about 220 yr, and 18.6 from 22 yr about 120 yr. The instrumental record
therefore cannot decide between these candidate periods. Two hypotheses make different predictions that long
records can test. An externally paced signal keeps a stable period and, at every site, a constant phase relative
to its clock (the phase may differ between sites). Internal (red-noise) variability shows a period that wanders
and a phase that drifts.

### B.2 Data

All proxy records were taken from the NOAA NCEI Paleoclimatology archive (plain-text NOAA templates). Six records
from independent sites were used for the test (Table B1). Ten further records that reuse the same corals or are
ENSO composites (Hon Tre wet-season SST, Fiji core 1F δ18O and Sr/Ca, Rarotonga δ18O, the Fiji–Tonga interdecadal
index, two Niño 3.4 reconstructions and a Niño cold-tongue reconstruction, and the Ocean2k western and eastern
Pacific coral composites) were analysed but never counted as independent. Central England Temperature
(1659–2025), the PAGES2k global-mean temperature median (1000–2000) and ERSSTv6 (1850–2025) served as comparison
series. Annual means were formed in each record's own year convention. δ18O and raw Sr/Ca ratios (not the records
already calibrated to SST) were multiplied by −1 so that positive means warm. The δ18O series also carry a
salinity signal.

**Table B1.** Independent coral records.

| Site | Region | Proxy | Years | N (yr) | Source |
|---|---|---|---|---|---|
| Hon Tre, Vietnam | South China Sea | Sr/Ca SST | 1600–2010 | 411 | Goodkin et al. (2021) |
| Fiji (core AB) | SW Pacific | δ18O (temperature and salinity) | 1617–2001 | 385 | Linsley et al. (2006) |
| Amédée Island, New Caledonia | SW Pacific | Sr/Ca anomaly | 1649–1999 | 351 | DeLong et al. (2012, 2013) |
| Rarotonga, Cook Islands | S Pacific | Sr/Ca | 1727–1996 | 270 | Linsley et al. (2000, 2006); Ren et al. (2003) |
| Dry Tortugas, Florida | Gulf of Mexico (Atlantic) | Sr/Ca SST | 1734–2008 | 275 | DeLong et al. (2014) |
| Northern Ryukyu Islands, Japan | NW Pacific | Sr/Ca SST | 1579–2008 | 430 | Kawakubo et al. (2017) |

### B.3 Methods

The main tests were specified before the data were analysed and were applied identically to every record and to
every surrogate series; the surrogate-based and common-period versions of the cross-site test were added as
sensitivity checks.

1. *Pre-treatment.* A LOESS trend (degree 2, span 100 yr) was removed to keep periods shorter than about 50 yr;
   the residual was standardised.
2. *Spectral test.* Periodogram with a split-cosine-bell taper (10 % at each end), zero-padded to four times the
   record length. The largest peak in the 15–25-yr band was compared with 2000 AR(1) surrogates whose lag-1
   autocorrelation and variance were estimated from the record before detrending; each surrogate was passed through
   the same pre-treatment. The decision is based on the band test (observed band maximum against surrogate band
   maxima), which accounts for the peak having been selected from the band. Because the undetrended records are
   redder (lag-1 autocorrelation 0.57–0.74 at the six sites, 0.34–0.45 after detrending), this null is
   conservative; with AR(1) parameters estimated after detrending, an anti-conservative choice, the six sites give
   band p = 0.14–1.0.
3. *Stationarity.* The record was band-passed (15–25 yr) and analysed in sliding 100-yr windows (step 10 yr). In
   each window the best-fitting period (grid 15–25 yr, step 0.1 yr) and amplitude were estimated. The spread of
   window periods was compared with 1000 AR(1) surrogates.
4. *Phase locking.* In each window the least-squares phase was estimated relative to three clocks: (a) the lunar
   nodal cycle, 18.613 yr, phase reference 2006.0 (the nodal modulation of high tides peaked in 2006 where diurnal
   tides dominate, Haigh et al. 2011; the major lunar standstill was in mid-2006, a constant offset that does not
   affect the tests); (b) a fixed 22.0-yr sinusoid; (c) the Hale cycle, represented by the annual mean SILSO
   sunspot number with its sign reversed at every cycle minimum (minima located on the 13-month running mean;
   1749–2025), band-passed like the records; its phase was estimated by regression on this series and its Hilbert
   quadrature. Because the sign assignment is arbitrary, only the constancy and the cross-site agreement of the
   Hale phase are interpreted. The constancy of the window phases (mean resultant length R) was tested against the
   AR(1) surrogates. In addition, the maximum correlation with the Hale series at lags of 0–5 yr was tested.
5. *Agreement across sites.* The full-record phases of the six independent sites were tested for clustering with
   a Rayleigh test (Mardia and Jupp 2000), using both the analytic approximation and surrogate phases, and for the
   common period 1750–1996.
6. *Dating uncertainty.* Each record was distorted by a random-walk age error that is zero at the youngest end and
   has a step standard deviation of 0.1 yr per year (standard deviation 1 yr after 100 yr, 2 yr after 400 yr). The
   window and full-record phase statistics were recomputed for 200 such realisations and compared with the 95 %
   surrogate thresholds of the undistorted analysis; the spectral test was not repeated.

### B.4 Results

- **No significant bidecadal peak.** None of the 19 series has a band-test p below 0.05 (range 0.14–1.0). For
  the six independent sites p = 0.46–1.0. Their peak periods (21.1, 22.0, 19.0, 24.0, 17.5 and 17.9 yr) scatter
  over the band instead of clustering at one value. Pointwise percentiles at the selected peak reach ≥ 95 % in six
  of the 19 series (among the independent sites only Hon Tre, 96.5 %), but they are biased high because the peak
  is selected from the band.
- **No stable period.** In 100-yr windows the best period (searched between 15 and 25 yr) moves across most of
  this range in most records and often reaches its edges. At no site is the period steadier than in red noise (New
  Caledonia comes closest, p = 0.09); at Rarotonga and Dry Tortugas it wanders even more than in red noise (spread
  larger than in 99 % and 93 % of the surrogates).
- **No phase lock at the individual sites.** Of 76 per-record tests (19 series × 4 statistics), 3 gave p < 0.05
  (a second Fiji core against the nodal cycle, p = 0.021; ERSSTv6 against the Hale series, p = 0.030; PAGES2k,
  maximum lagged correlation with the Hale series, p = 0.047), close to the 3.8 expected by chance. None of the six
  independent records is among them.
- **No common phase across sites.** Rayleigh p = 0.13 for the nodal cycle, 0.38 for the 22-yr sinusoid and 0.21
  for the Hale cycle (surrogate-based 0.12, 0.38, 0.21; common period 1750–1996: 0.12, 0.46, 0.22). The nodal
  phases cluster most (five of six within −5° to +83°), but not significantly. Pooling all 16 coral and proxy
  records gives p = 0.006 for the Hale phase, but these records share corals and ENSO teleconnections and are not
  independent; we do not regard this as evidence.
- **Dating errors of the assumed size** change the window-phase statistic R by less than 0.05 in every record, but
  full-record phases, which enter the cross-site test, are less robust (circular spread up to about ±50° for the
  longest record, Hon Tre).

**Table B2.** Results for the six independent sites.

| Site | Peak (yr) | Band p | Window period, median (s.d.) | p_stable | p nodal | p 22-yr | p Hale | p lag-corr. | Full-record phase nodal / 22-yr / Hale (°) |
|---|---|---|---|---|---|---|---|---|---|
| Hon Tre | 21.1 | 0.46 | 20.8 (2.9) | 0.58 | 0.72 | 0.89 | 0.85 | 0.77 | −5 / 164 / 78 |
| Fiji AB | 22.0 | 1.00 | 23.7 (2.8) | 0.51 | 0.53 | 0.49 | 0.07 | 0.11 | 24 / −5 / −55 |
| New Caledonia | 19.0 | 0.54 | 18.4 (1.6) | 0.09 | 0.20 | 0.26 | 0.12 | 0.52 | 64 / 39 / 22 |
| Rarotonga | 24.0 | 0.96 | 16.9 (4.5) | 0.99 | 0.49 | 0.65 | 0.91 | 0.46 | 63 / −94 / −131 |
| Dry Tortugas | 17.5 | 0.99 | 17.2 (3.9) | 0.93 | 0.74 | 0.64 | 0.94 | 0.97 | 83 / −57 / −20 |
| Ryukyu | 17.9 | 0.67 | 18.0 (2.4) | 0.30 | 0.60 | 0.79 | 0.36 | 0.24 | −96 / −82 / −22 |
| Rayleigh p, 6 sites (surrogate) | | | | | 0.13 (0.12) | 0.38 (0.38) | 0.21 (0.21) | | |

p_stable = fraction of AR(1) surrogates with a smaller spread of window periods (small = steadier than noise);
phases relative to 2006.0 (nodal), 2000.0 (22-yr) and the Hale-signed series (sign arbitrary).

### B.5 Limitations

- Six sites are few, and some share regional climate. A weak lock with the amplitude seen in the global mean
  (about 0.03 K) would be hard to detect in noisy local proxies: replicate Sr/Ca records from one colony differ by
  about 0.4 °C for monthly values (DeLong et al. 2013), and interannual SST anomalies at the Ryukyu site are only
  about 0.5 °C (Kawakubo et al. 2017). Failure to reject therefore does not prove absence.
- The assumed age error is optimistic for single-core chronologies: replicated cores from New Caledonia revealed
  chronological discrepancies of −2.3 and +3.7 yr per century in single-core reconstructions (DeLong et al. 2013).
  Errors of 3–5 yr in the older sections would blur a phase lock at 18.6–22 yr beyond detection. The New Caledonia
  and Dry Tortugas chronologies are replicated and checked by U–Th dating; the Fiji and Rarotonga chronologies rest
  on counting seasonal cycles in single cores, and the data files for Hon Tre and the Ryukyus give no dating
  uncertainty.
- A common clock need not produce the same phase at all sites: the local SST response can differ in sign and lag
  (the bidecadal pattern in ERSSTv6 has opposite signs in the tropical Pacific and in the central North and
  subtropical south-west Pacific). The cross-site Rayleigh test therefore has power only for an in-phase response;
  the per-site window tests do not require it.
- 100-yr windows cannot separate 18.6, 20 and 22 yr. They test whether the phase stays constant, not which clock
  applies.
- The nodal tide peaks at different dates in diurnal (2006, 2024) and semi-diurnal (1997, 2015) tidal regimes
  (Haigh et al. 2011). The cross-site test assumed one phase everywhere and could miss a regime-dependent lock.
- The Hale series begins in 1749, so the Hale tests use only part of each record.
- Only AR(1) surrogates were used; long-memory noise would make the results even less significant.

### B.6 Conclusion

Long coral records from six independent sites show no stationary, significant ~20-yr signal and no detectable
consistent phase relation to the lunar nodal cycle, the Hale cycle or a fixed 22-yr period. The behaviour is
consistent with internal red-noise variability. A weak external pacing, or a lock whose phase differs between
sites or tidal regimes, cannot be excluded with the present data.

### B.7 References

- Clette, F., Lefèvre, L. (2016). The new sunspot number: assembling all corrections. *Solar Physics* 291(9–10),
  2629–2651. doi:10.1007/s11207-016-1014-y
- DeLong, K. L., Quinn, T. M., Taylor, F. W., Lin, K., Shen, C.-C. (2012). Sea surface temperature variability in
  the southwest tropical Pacific since AD 1649. *Nature Climate Change* 2(11), 799–804. doi:10.1038/nclimate1583
- DeLong, K. L., Quinn, T. M., Taylor, F. W., Shen, C.-C., Lin, K. (2013). Improving coral-based paleoclimate
  reconstructions by replicating 350 years of coral Sr/Ca variations. *Palaeogeography, Palaeoclimatology,
  Palaeoecology* 373, 6–24. doi:10.1016/j.palaeo.2012.08.019
- DeLong, K. L., Flannery, J. A., Poore, R. Z., Quinn, T. M., Maupin, C. R., Lin, K., Shen, C.-C. (2014). A
  reconstruction of sea surface temperature variability in the southeastern Gulf of Mexico from 1734 to 2008 C.E.
  using cross-dated Sr/Ca records from the coral *Siderastrea siderea*. *Paleoceanography* 29(5), 403–422.
  doi:10.1002/2013PA002524
- Freund, M. B., Henley, B. J., Karoly, D. J., McGregor, H. V., Abram, N. J., Dommenget, D. (2019). Higher
  frequency of Central Pacific El Niño events in recent decades relative to past centuries. *Nature Geoscience*
  12(6), 450–455. doi:10.1038/s41561-019-0353-3
- Goodkin, N. F., Samanta, D., Bolton, A., Ong, M. R. K., Hoang, P. K., Vo, S. T., Karnauskas, K. B., Hughen, K. A.
  (2021). Natural and anthropogenic forcing of multi-decadal to centennial scale variability of sea surface
  temperature in the South China Sea. *Paleoceanography and Paleoclimatology* 36(10), e2021PA004233.
  doi:10.1029/2021PA004233
- Haigh, I. D., Eliot, M., Pattiaratchi, C. (2011). Global influences of the 18.61 year nodal cycle and 8.85 year
  cycle of lunar perigee on high tidal levels. *Journal of Geophysical Research* 116(C6), C06025.
  doi:10.1029/2010JC006645
- Huang, B., et al. (2025). Extended Reconstructed Sea Surface Temperature version 6 (ERSSTv6), Parts I and II.
  *Journal of Climate* 38(4), 1105–1121 and 1123–1136. doi:10.1175/JCLI-D-23-0707.1 (Part I)
- Kawakubo, Y., Alibert, C., Yokoyama, Y. (2017). A reconstruction of subtropical western North Pacific SST
  variability back to 1578, based on a *Porites* coral Sr/Ca record from the northern Ryukyus, Japan.
  *Paleoceanography* 32(12), 1352–1370. doi:10.1002/2017PA003203
- Linsley, B. K., Wellington, G. M., Schrag, D. P. (2000). Decadal sea surface temperature variability in the
  subtropical South Pacific from 1726 to 1997 A.D. *Science* 290(5494), 1145–1148.
  doi:10.1126/science.290.5494.1145
- Linsley, B. K., Kaplan, A., Gouriou, Y., Salinger, J., deMenocal, P. B., Wellington, G. M., Howe, S. S. (2006).
  Tracking the extent of the South Pacific Convergence Zone since the early 1600s. *Geochemistry, Geophysics,
  Geosystems* 7(5), Q05003. doi:10.1029/2005GC001115
- Linsley, B. K., Zhang, P., Kaplan, A., Howe, S. S., Wellington, G. M. (2008). Interdecadal-decadal climate
  variability from multicoral oxygen isotope records in the South Pacific Convergence Zone region since 1650 A.D.
  *Paleoceanography* 23(2), PA2219. doi:10.1029/2007PA001539
- Manley, G. (1974). Central England temperatures: monthly means 1659 to 1973. *Quarterly Journal of the Royal
  Meteorological Society* 100(425), 389–405. doi:10.1002/qj.49710042511
- Mardia, K. V., Jupp, P. E. (2000). *Directional Statistics*. Wiley, Chichester. doi:10.1002/9780470316979
- PAGES 2k Consortium (2019). Consistent multidecadal variability in global temperature reconstructions and
  simulations over the Common Era. *Nature Geoscience* 12(8), 643–649. doi:10.1038/s41561-019-0400-0
- Parker, D. E., Legg, T. P., Folland, C. K. (1992). A new daily Central England Temperature series, 1772–1991.
  *International Journal of Climatology* 12(4), 317–342. doi:10.1002/joc.3370120402
- Ren, L., Linsley, B. K., Wellington, G. M., Schrag, D. P., Hoegh-Guldberg, O. (2003). Deconvolving the δ18O
  seawater component from subseasonal coral δ18O and Sr/Ca at Rarotonga in the southwestern subtropical Pacific
  for the period 1726 to 1997. *Geochimica et Cosmochimica Acta* 67(9), 1609–1621.
  doi:10.1016/S0016-7037(02)00917-1
- Tierney, J. E., Abram, N. J., Anchukaitis, K. J., Evans, M. N., Giry, C., Kilbourne, K. H., Saenger, C. P., Wu,
  H. C., Zinke, J. (2015). Tropical sea surface temperatures for the past four centuries reconstructed from coral
  archives. *Paleoceanography* 30(3), 226–252. doi:10.1002/2014PA002717
- Wilson, R., Cook, E., D'Arrigo, R., Riedwyl, N., Evans, M. N., Tudhope, A., Allan, R. (2010). Reconstructing
  ENSO: the influence of method, proxy data, climate forcing and teleconnections. *Journal of Quaternary Science*
  25(1), 62–78. doi:10.1002/jqs.1297

Data: NOAA NCEI Paleoclimatology, study IDs 32392, 16216, 13035, 6089, 16217, 34372 (independent sites) and 6216,
11749, 26270, 17955 (dependent records); HadCET (Manley 1974; Parker et al. 1992); PAGES 2k Consortium (2019);
ERSSTv6 (Huang et al. 2025). Sunspot numbers: WDC-SILSO, Royal Observatory of Belgium, Brussels (version 2.0;
Clette and Lefèvre 2016). Code: R scripts (R 4.6.2) in this folder. Figure B1: periodograms with the AR(1)
95 % level (`fig_spectra.png`); Figure B2: Morlet wavelet power, no significance contours computed (`fig_wavelet.png`); Figures B3 and B4: period, amplitude and phase in sliding 100-yr windows (`fig_sliding_period_amp.png`, `fig_sliding_phase.png`).
