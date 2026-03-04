# CosmoBurial

MATLAB tools for cosmogenic nuclide burial dating using paired ¹⁰Be and ²⁶Al concentrations. The suite includes three levels of model complexity to accommodate different burial depths, depositional histories, and available geological constraints.

---

## Dependencies

These functions require the **CRONUS-Earth cosmogenic nuclide calculators**, specifically:
- `LSD.m` — LSD scaling framework (Lifton et al. 2014)
- `fit_Pmu_with_exp_1A.m` — muon production fitting (Balco 2017)

Both are available at the [CRONUS-Earth calculator repository](https://bitbucket.org/cronusearth/cronus-calc).

---

## Functions

### `SimpleBurial_nomuons.m`
Computes burial age and pre-burial erosion rate assuming instantaneous burial. Post-burial cosmogenic production is ignored. Appropriate for burial depths greater than ~500 g/cm².

```matlab
[Tburial, Tburial_unc, Ero, Ero_unc] = SimpleBurial_nomuons(...
    be10, be10_unc, al26, al26_unc, lat, lon, elev)
```

---

### `SimpleBurial_muons.m`
Same as above but includes post-burial muon production at the burial depth as a single lumped term. Use when burial depth is shallow enough that the muon flux at depth is non-negligible.

```matlab
[Tburial, Tburial_unc, Ero, Ero_unc] = SimpleBurial_muons(...
    be10, be10_unc, al26, al26_unc, lat, lon, elev, mass_depth)
```

> **Note:** At burial depths greater than ~2000 g/cm², the lumped muon term can destabilize the solver due to the muon Al/Be ratio at depth exceeding the surface spallation ratio. Use `SimpleBurial_muons2exp` in those cases.

---

### `SimpleBurial_muons2exp.m`
Extends the simple burial model with a full two-exponential muon production scheme (Balco 2017). Inherited muon production and post-burial muon accumulation are each treated as two independent exponential components with distinct attenuation lengths. This avoids the ratio instability of the single-term approach and is the recommended muon-corrected model for most applications.

```matlab
[Tburial, Tburial_unc, Ero, Ero_unc] = SimpleBurial_muons2exp(...
    be10, be10_unc, al26, al26_unc, lat, lon, elev, mass_depth)
```

---

### `ComplexBurial_muons.m`
Full burial model that accounts for post-burial cosmogenic production during sediment accumulation. The burial depth increases with time as sediment accumulates at rate `Ar`, and both spallogenic and muogenic production are integrated analytically over the burial history. Requires `GetProductionRatesBurial.m` and an initial age estimate from `SimpleBurial_nomuons`.

```matlab
[Burial_age, Age_unc, PaleoEro, Ero_unc] = BurialCalcComplex(...
    be10, be10_unc, al26, al26_unc, lat, lon, ...
    sedrateflag, sedrate, sink_elev, mass_depth, density)
```

`sedrateflag = 1`: provide a fixed sedimentation rate via `sedrate` [cm/yr].  
`sedrateflag = 0`: sedimentation rate is derived internally from the simple burial age estimate and the known mass depth.

---

## References

Balco, G. (2017). Production rate calculations for cosmic-ray muon-produced ¹⁰Be and ²⁶Al benchmarked against geological calibration data. *Quaternary Geochronology*, 39, 150–173.

Basunia, M.S. & Hurst, A.M. (2016). Nuclear Data Sheets for A = 26. *Nuclear Data Sheets*, 134, 1–148.

Borchers, B. et al. (2016). Geological calibration of spallation production rates in the CRONUS-Earth project. *Quaternary Geochronology*, 31, 188–198.

Chmeleff, J. et al. (2010). Determination of the ¹⁰Be half-life by multicollector ICP-MS and liquid scintillation counting. *Nuclear Instruments and Methods in Physics Research B*, 268, 192–199.

Granger, D.E. & Muzikar, P.F. (2001). Dating sediment burial with in situ–produced cosmogenic nuclides: theory, techniques, and limitations. *Earth and Planetary Science Letters*, 188, 269–281.

Lifton, N. et al. (2014). Scaling in situ cosmogenic nuclide production rates using analytical approximations to atmospheric cosmic-ray fluxes. *Earth and Planetary Science Letters*, 386, 149–160.

---

## Authors

Helbert Garcia-Delgado & Benjamin Guerrero  
Syracuse University — Department of Earth and Environmental Sciences
