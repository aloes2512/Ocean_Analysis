# The coral experiment: is the ~20-year line in the ocean temperature paced by the sun or the moon?

This folder documents one experiment done on 22 November 2026 during the re-analysis of the Ocean_Analysis paper. It
contains a plain-language explanation (this file), a draft section for the manuscript (`manuscript_section.md`), the
data, the R scripts, the result tables and the figures. Everything can be re-run (see "How to re-run" below).

## 1. The question

The global mean sea surface temperature (ERSSTv6, 1850–2025) contains a weak oscillation with a period of about
20 years and an amplitude of about 0.03 K. Three candidate explanations have very similar periods:

- the **lunar nodal cycle** (18.6 years), which modulates the tides;
- the **Hale cycle** of the sun (about 22 years, two sunspot cycles with reversing magnetic polarity);
- **internal variability** of the ocean (for example the Pacific decadal oscillation), which needs no external
  clock at all and is excited by random weather.

The thermometer record is 176 years long. That is too short to tell 18.6, 20 and 22 years apart: two periods can
only be separated if the record is longer than 1/|1/T1 − 1/T2|. That gives about 270 years for 18.6 vs 20 years,
220 years for 20 vs 22 years and 120 years for 18.6 vs 22 years, and somewhat more when there is noise.

## 2. The idea of the test

The three explanations make different predictions over several centuries:

- If an **external clock** (moon or sun) drives the line, it should keep **the same period** for centuries, and at
  each place it should keep **a fixed phase** relative to that clock. The phase may differ from place to place,
  because the local response can have a different sign or delay, but it must not drift.
- If the line is **internal variability**, its period and its phase **wander** from century to century, like
  filtered random noise.

So the test needs temperature records that are much longer than 176 years. Thermometers do not go back that far in
the ocean. Corals do.

## 3. Why corals, and which ones

Reef corals build their skeleton (aragonite) in yearly density bands that are visible in X-ray images, similar to
tree rings. While the skeleton grows, its chemistry records the water temperature:

- the **strontium-to-calcium ratio (Sr/Ca)** is lower in warmer water;
- the **oxygen isotope ratio (δ18O)** is lower in warmer or fresher water, so it also carries a salinity signal.

Measured in small steps along the growth axis and dated by counting the yearly cycles, this gives a temperature
series for the place where the coral lives. Large colonies, mostly of the genus *Porites*, live for several hundred
years. All records used here are public in the NOAA NCEI Paleoclimatology archive.

We used **six records from six independent sites**, each 270–430 years long:

| Site | Region | Proxy | Years | Length | Source |
|---|---|---|---|---|---|
| Hon Tre, Vietnam | South China Sea | Sr/Ca SST | 1600–2010 | 411 yr | Goodkin et al. 2021 |
| Fiji (core AB) | SW Pacific | δ18O | 1617–2001 | 385 yr | Linsley et al. 2006 |
| Amédée Island, New Caledonia | SW Pacific | Sr/Ca | 1649–1999 | 351 yr | DeLong et al. 2012, 2013 |
| Rarotonga, Cook Islands | S Pacific | Sr/Ca | 1727–1996 | 270 yr | Linsley et al. 2000, 2006 |
| Dry Tortugas, Florida | Gulf of Mexico | Sr/Ca SST | 1734–2008 | 275 yr | DeLong et al. 2014 |
| Northern Ryukyu Islands, Japan | NW Pacific | Sr/Ca SST | 1579–2008 | 430 yr | Kawakubo et al. 2017 |

Ten further records that reuse the same corals or combine several of them (for example El Niño reconstructions)
were analysed too but not counted as independent evidence. Central England Temperature (since 1659), the PAGES2k
global temperature reconstruction and ERSSTv6 were run through the same steps for comparison.

Note that a coral measures the temperature at one place, not the global mean. For the question "is there a fixed
phase to an external clock?" this is not a problem: a global clock would have to show up with a constant phase at
each place.

## 4. What exactly was done (the same steps for every record)

1. **Remove the slow trend.** A smooth curve (LOESS, span 100 years) was subtracted, so that only periods shorter
   than about 50 years remain. The rest was scaled to unit variance.
2. **Is there a significant peak between 15 and 25 years?** We computed the spectrum (periodogram) and compared the
   highest peak in the 15–25-year band with 2000 artificial red-noise series (AR(1)) that have the same
   year-to-year memory and variance as the record. The test takes into account that we picked the highest peak out
   of the band.
3. **Is the period stable?** We cut the record into 100-year windows (shifted by 10 years each time) and fitted the
   best period between 15 and 25 years in each window. If the spread of these periods is smaller than in red noise,
   the period is steadier than chance.
4. **Is there a fixed phase to a clock?** In each window we measured the phase of the oscillation relative to three
   clocks: the lunar nodal cycle (18.613 years), a fixed 22-year sine wave, and the Hale cycle (the yearly sunspot
   number with its sign flipped at every solar minimum). If the phase stays constant from window to window more
   than in red noise, that is a sign of phase locking.
5. **Do the sites agree?** We tested whether the full-record phases of the six independent sites cluster together
   (Rayleigh test).
6. **What about dating errors?** Old coral sections can be misdated by a few years. We distorted each record 200
   times with random dating errors and repeated the phase tests.

All tests were fixed before the data were analysed and were applied the same way to every record.

## 5. What came out

- **No significant ~20-year peak** at any of the six sites (p = 0.46–1.0; a significant result would need p
  below 0.05). The peak periods (21.1, 22.0, 19.0, 24.0, 17.5 and 17.9 years) scatter over the whole band instead
  of meeting at one value.
- **The period is not stable.** Within one record it moves through most of the 15–25-year range from century to
  century. At no site is it steadier than red noise. At Rarotonga and Dry Tortugas it wanders even more than red
  noise.
- **No fixed phase** at any independent site, neither to the moon, nor to the 22-year sine, nor to the Hale cycle.
  Of all 76 single tests (19 series × 4 tests), 3 gave p below 0.05. About 4 are expected by pure chance, and none
  of the three is an independent site.
- **The sites do not share a phase** (p = 0.13 for the moon, 0.38 for 22 years, 0.21 for Hale).

**Conclusion.** Over the last 300–400 years the ~20-year oscillation in these corals behaves like internal ocean
variability, not like a signal paced by the sun or the moon. This does not prove that there is no solar or lunar
influence: a weak one could be hidden in the noise of local corals (the global signal is only 0.03 K), and the
dating of single cores can be off by 2–4 years per century, which would blur a lock at these periods. The limits
are listed in full in `manuscript_section.md`, section B.5.

## 6. How this relates to Alf's current goal

Alf wrote that these analyses go well beyond his goal, which is to extract similar dominant quasi-harmonics from
different regions (near the equator such as ENSO and the Indian Ocean; middle latitudes such as the Mediterranean,
the Gulf of Mexico and the Pacific; later perhaps northern latitudes) and then derive coupling factors with dynamic
indices such as El Niño. This folder is not meant to change that direction. It documents what was done on
22 November, because the earlier emails did not explain it well.

Two points may still be useful for the regional work:

- Several of these corals sit in the regions Alf named: the Gulf of Mexico (Dry Tortugas), the tropical and
  south-west Pacific (Fiji, New Caledonia, Rarotonga) and the western Pacific (South China Sea, Ryukyu). They are
  regional temperature records that reach back 270–430 years, so they can show whether a quasi-harmonic found in
  the 176-year thermometer data keeps its period over a longer time.
- The same tests (red-noise comparison, stability in windows) can be applied to any regional series Alf extracts.

## 7. Files

| File | Content |
|---|---|
| `README.md` | this explanation |
| `manuscript_section.md` | the same experiment written as a manuscript section (data, methods, results, limitations, references) |
| `fig_spectra.png` | spectra of all records with the red-noise 95 % level (Figure B1) |
| `fig_wavelet.png` | wavelet power: which periods are strong at which time (Figure B2, no significance contours) |
| `fig_sliding_period_amp.png` | best period and amplitude in sliding 100-year windows (Figure B3) |
| `fig_sliding_phase.png` | phase relative to the three clocks in sliding 100-year windows (Figure B4) |
| `out/table_records.csv` | all test results per record |
| `out/table_fullphases.csv` | full-record phases (degrees) used in the cross-site test |
| `out/table_dating.csv` | effect of the simulated dating errors |
| `out/run_log.txt`, `out/post_log.txt` | console output of the last run |
| `A2_load.R`, `A2_lib.R`, `A2_run.R`, `A2_post.R` | the R scripts |
| `data/` | the input data (NOAA NCEI coral files with their metadata in `data/meta/`, HadCET, PAGES2k, SILSO sunspots, the global ERSSTv6 series) |

## 8. How to re-run

Requires R (tested with 4.6.2) with the base packages only (`parallel` is part of base R). In R:

```r
setwd("path/to/Ocean_Analysis/coral_test")
source("A2_load.R")   # reads the data          (seconds)
source("A2_run.R")    # all tests, 2000 surrogates per record (about 4 minutes on 7 cores)
source("A2_post.R")   # tables and figures      (seconds)
```

The random seeds are fixed, so the numbers are reproducible.

## Data credits

Coral records: NOAA NCEI Paleoclimatology (study IDs 32392, 16216, 13035, 6089, 16217, 34372; dependent records
6216, 11749, 26270, 17955), see the references in `manuscript_section.md`. Central England Temperature: Met Office
Hadley Centre (HadCET). PAGES 2k Consortium (2019). ERSSTv6: NOAA (Huang et al. 2025). Sunspot numbers: WDC-SILSO,
Royal Observatory of Belgium, Brussels.
