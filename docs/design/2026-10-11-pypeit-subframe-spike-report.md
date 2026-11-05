---
title: "🔭 mastpipe — PypeIt per-fiber spike: results"
subtitle: "Throwaway experiment · 2026-10-11"
---

<div class="lede">
Follow-up to the question "why not give PypeIt one sub-frame per trace?". A throwaway PypeIt spectrograph class was written that serves every fiber as its own detector, and run on the three-fiber laboratory Th-Ar data of the HighSpec laboratory reduction notebook and on a simulated 45-trace survey frame. Spike by Claude 🤖 at the project lead's request; numbers below are measured, not estimated. Input for the team meeting.
</div>

## ✅ Short answer

**It works, with caveats.** PypeIt 2.0.1 can take each fiber as a bounded-detector long slit, solve the wavelength scale per fiber over the full 14 nm from a Th-Ar template, and extract spectra with correct errors, for 3 lab fibers and for 45 simulated fibers, with no failures once about fifteen parameters are set correctly. What it does not provide is unchanged: sky from sky fibers, the arc-fiber drift transfer, combination, and the observatory service. Its wavelength residuals are 4 to 10 times the notebook's physical-model solution, and configuring it for fibers took seven debugging rounds against defaults made for slits.

## 🧪 What was built

| Piece | Detail |
|---|---|
| Environment | fresh conda env, PypeIt **2.0.1**, Python 3.12 |
| Spectrograph class | ~150 lines; `ndet` = number of fibers; `get_rawimage` returns a horizontal sub-frame around each trace; layout from a JSON environment variable |
| Wavelength template | UVES Th-Ar atlas (Murphy et al.) convolved to R = 23 000, mapped onto the pixel grid with the notebook's grating constants, converted to vacuum; dispersion direction found by cross-correlation with the lab arc |
| Lab data | The notebook's bias, flat, Th-Ar focus frame as arc, a second focus frame as "science" |
| Synthetic frame | 45 Gaussian traces at 8.5 px pitch, 3 px FWHM, 5 arc / 20 science / 20 sky fibers, flat, arc, science with absorption lines, truth saved |

## 📊 Results

**Laboratory data, 3 fibers (sub-frames of ±9 rows):**

| Quantity | PypeIt | Laboratory notebook |
|---|---|---|
| Wavelength solution | 48 / 56 / 48 lines, rms **0.059 / 0.112 / 0.113 px** (0.5 / 1.0 / 1.0 pm) | 23 / 25 / 19 lines, rms 0.016 / 0.038 / 0.010 px (0.14 / 0.33 / 0.09 pm) |
| Solved coverage, fiber 0 (air) | 513.499 – 527.310 nm | 513.498 – 527.292 nm |
| Fitted spatial FWHM | 2.9 – 3.1 px | 2.58 px (line width along dispersion) |
| Extraction | boxcar consistent with a direct 19-row sum; optimal extraction S/N 50 – 60 per pixel on the arc | aperture sum |
| Run time | **25 s** for the whole chain, 3 fibers; **≈4 s** per extra science frame | about a minute for the full notebook |

**Simulated survey frame, 45 fibers:**

| Quantity | Sub-frame ±6 rows | Sub-frame ±4 rows |
|---|---|---|
| Wavelength solutions | 45 of 45; rms median **0.119 px**, range 0.078 – 0.179; 49 – 61 lines | same |
| Extracted / true flux, science fibers | median **0.93** (0.64 – 1.00) | 0.71 (apertures truncated) |
| Error bars | residual scatter 22.7 vs reported 21.8 counts: **accurate** | — |
| Crosstalk indicator (sky fibers vs science shape) | correlation 0.07 | 0.08 |
| Run time, full night (3 bias, flat, arc, 1 science) | **314 s** | 286 s |

At 1 pm per pixel, PypeIt's wavelength rms corresponds to about 0.6 km s⁻¹ at 518 nm; the notebook's to 0.05 – 0.2 km s⁻¹. Both meet a few-km s⁻¹ goal; the physical model is clearly the better calibration and could still be applied to PypeIt's extracted arcs.

## 🧱 What had to be overridden

Each of these cost one failed run before it was found:

- 🔧 `use_overscan = False` — the camera has no overscan; PypeIt assumes one.
- 🔧 `slitedges.edge_thresh = 1e6` — PypeIt's edge detector **traced the fiber image itself as a 4-pixel slit**, truncating every aperture; edge detection must be switched off so the bounded-detector fallback is used.
- 🔧 `slitedges.det_buffer = 1` — the default buffer of 5 px makes the bounded-detector code crash on narrow sub-frames.
- 🔧 `reduce.trim_edge = 0,0` and `findobj.find_trim_edge = 0,0` — default trims of 3 and 5 px per side would eat a fiber sub-frame.
- 🔧 `findobj.skip_skysub = True`, `skysub.no_local_sky = True` — no sky model is meaningful inside a fiber.
- 🔧 Manual extraction per detector (`det:spat:spec:fwhm:boxcar_radius`) — the automatic object finder found nothing on the emission-line "object"; positions come from the layout table, which is what the design proposes anyway.
- 🔧 `wavelengths.refframe = observed` — a heliocentric correction was applied by default to laboratory data.
- 🔧 Template file format: one row per pixel, a `BINSPEC` header card, vacuum wavelengths; none of this is documented for user templates.
- 🔧 Flat-field use disabled (`use_pixelflat`, `use_illumflat`) — the illumination model would divide out the fiber profile.

None of these is a bug; they are slit assumptions that fibers violate. The point is that every one is silent until the output is checked against a direct sum or a truth.

## ⚖️ What this changes in the comparison

| | Independent package (current plan) | PypeIt engine per fiber (this spike) |
|---|---|---|
| Calibration, traces, extraction | port the notebook's code; write Horne extraction (~200 lines) | provided, mature, with cosmic-ray rejection and accurate errors |
| Wavelength solution | The notebook's grating model, 0.1 – 0.3 pm | generic polynomial on PypeIt's line list, ≈1 pm; the grating model could be run afterwards on the extracted arcs |
| Sky from sky fibers, arc-fiber transfer, combination, products, QA pages, service | ours | **still ours**, operating on PypeIt's `spec1d` outputs |
| Header changes | one YAML file | the spectrograph class |
| Quick-look run model | our runner, idempotent, one failure = one frame | `pypeit_setup` + `run_pypeit`, ~5 min of calibrations per night plus ~1 min per exposure for 45 fibers, many intermediate files |
| Dependencies | numpy, scipy, astropy, pandas, matplotlib, pydantic, pyyaml | adds PypeIt with PyQt6, ginga, scikit-learn, bottleneck, …; Windows "most installation issues" per PypeIt docs, native or WSL |
| Risk | we own every line | fiber behaviour depends on defaults made for slits; each PypeIt release needs re-validation |
| Plan 1 effort | as planned | shrinks: skeleton, ingest, layout, simulator, products stay; calibration/trace/extraction tasks go |

## 🧭 Recommendation for the meeting

Present both as viable, with these numbers, and decide on priorities rather than feasibility:

1. If the group values **a maintained extraction engine and less code to own**, choose "PypeIt per fiber inside mastpipe": the spike class becomes the adapter, and mastpipe keeps ingest, layout, drift transfer, sky, combination, products and the service.
2. If the group values **a lean, fully controlled service on the observatory machine and the better wavelength model**, keep the independent plan and borrow PypeIt's Horne extraction code only.

Claude's lean, for the record: option 2 for the observatory quick-look, with the spike's PypeIt setup kept as an independent cross-check during commissioning. The deciding factors were the seven silent slit assumptions and the dependency weight, not capability.

## 📁 Where things are

Spike scripts, template, run file, result tables and a README with every setting are saved under the project folder, `spike/2026-10-11-pypeit-subframes/` (git-ignored, throwaway). The conda environment `pypeit2` holds the patched PypeIt install.
