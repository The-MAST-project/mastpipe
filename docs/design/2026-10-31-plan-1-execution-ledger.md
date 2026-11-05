# SDD ledger — plan: docs/design/2026-10-29-mastpipe-plan-1-highspec.md

Spec: docs/design/2026-09-24-mastpipe-design.md + docs/design/2026-10-29-mastpipe-design-revision-1.md (revision wins).
Repo: fresh local git repo on `main` (plan Task 1 mandates `git init -b main`, local only, no remote). Base commit = .gitignore only.
Ruling: no worktree/feature branch — the repo is new and the plan commits to `main` locally; the user approved the plan — costs nothing if wrong (no shared branch exists).
Ruling: every subagent runs on the session model (fable) — repo CLAUDE.md "One model throughout" and memory feedback-mastpipe-workflow override the skill's cheap-model guidance — costs tokens, not correctness.

## Pre-flight scan

| Pair / task | Produces vs consumes | Finding / ruling |
|---|---|---|
| T2 Spectrum.latest(upto) ↔ T16 combine(version=None) | upto caps the version; None = no cap = most advanced | consistent |
| T4 make_line_table.py ↔ T11 atlas functions | script imports T11; CSV needed by T4/T6 before T11 exists | plan: T4 ships CSV from a throwaway run; T11 regenerates and diffs. Ruling: T4 implementer may generate the CSV with a throwaway script run against the local atlas (path given only in the dispatch, never in a tracked file) — cost if wrong: a CSV diff at T11 |
| T6 simulator ↔ T4 load_line_table | line table columns wavelength_nm, species, logI, usable | consistent |
| T8 classify_arc / comb_shift ↔ T14, T18 | signatures reused by name | consistent |
| T11 model_from_solution ↔ T12/14/17/18 | rebuilds model from WavelengthSolution(kind, params, n_pix) + geom | consistent after self-review fix |
| T12 solve_full_arc(start=model) "slid against fiber 0" ↔ T13 anchor_trace param | anchor model belongs to anchor_trace, not necessarily fiber 0 | Ruling: add `ref_fiber: int = 0` to solve_full_arc (the fiber the single start model belongs to); T13 passes anchor_trace — cost if wrong: one extra parameter |
| T17 write_spectra(path, es, prov) "DRIFT table when drift arcs used" ↔ T18 runner | no argument carries the drift table | Ruling: `write_spectra(path, es, prov, drift: pd.DataFrame | None = None)`; runner passes the weighted per-fiber table (columns as DriftModel.per_fiber plus `weight`) — cost: one optional arg |
| T17 read_reference(path, geometry_lookup) ↔ T18 find_reference | library must supply geometry | Ruling: CalibLibrary resolves geometry through the registry (load_instrument → configuration geometry + camera pixel size → grating_geometry) — cost: library imports registry (pipeline layer, allowed) |
| T3 WavecalParams.min_arc_fibers_lit | consumed by nothing | Ruling: omit it (YAGNI) — cost: none |
| T9 fit_continuum raw flux untouched ↔ T12 shift_to_reference(flux) / T13 panels | callers use flux − continuum | consistent (stated in T9) |
| T5 arc_kind None ↔ T8 classify_arc ↔ T18 step 4 | runner fills arc_kind when None | consistent |
| T14 pair_arcs ↔ T18 options.time_policy | option overrides mode default | consistent |
| Review Focus 1–5 ↔ named tests | all five tests named in T18, T15, T14, T10, T5/T18 | consistent |
| Each task's own text | tests vs code vs files | T1–T19 self-consistent after plan self-review |

## Progress
Task 1: dispatched, BASE dd1f8d6 (implementer agent af44afa708f2ed827)
Task 1: complete (commits dd1f8d6..3ac58d2, review clean)
Task 1: minor (deferred): pyproject license table form is pre-PEP 639 (setuptools>=77 deprecation) — fix: license = "BSD-3-Clause", license-files, drop classifier, requires setuptools>=77
Task 1: minor (deferred): cli.build_parser discards the subparsers action — later tasks need a registration point
Task 1: minor (deferred): no automated test for scripts/privacy_scan.py
Task 1: minor (deferred): README "design documents live outside the repository" wording (they are git-ignored in docs/design/)
Task 1: minor (deferred): CLAUDE.md (draft copy) lints `src tests` not `src tests scripts` and mentions a pre-commit scan that is not installed — needs user approval to edit CLAUDE.md
Task 1: minor (deferred): tests/ packages lack __init__.py; duplicate basenames across tests/unit, e2e, regression would break collection
Ruling: fold Task 1 minors 1, 2 and 6 (PEP 639 metadata, subparser registration point, tests __init__.py) into the Task 2 dispatch as small additions — they are mechanical and the parser one is needed by Tasks 13/18 — costs a slightly larger Task 2 diff
Task 2: dispatched, BASE 3ac58d2
Task 2: complete (commits 3ac58d2..251ae82, review clean; ⚠️ COMBINED_VERSIONS=("raw",*VERSIONS) resolved: matches Task 16 plan text)
Task 2: minor (deferred): make_identifier does not normalise a tz offset to UTC (naive/Z fine)
Task 2: minor (deferred): TraceSet.by_telescope keys None for role=None and overwrites duplicate roles silently
Task 2: minor (deferred): cli.main args.func unguarded; CommandRegistrar aliases argparse._SubParsersAction at import time
Task 2: minor (deferred): test helpers placement; _check_length TypeError on 0-d arrays; raw-arrays-never-replaced is convention not enforcement (docstring note wanted)
Task 3: dispatched, BASE 251ae82
Task 3: complete (commits 251ae82..b6ba651, review clean)
Ruling: config_hash excludes EffectiveConfig.layers (hash over values only, so a restated override does not change the idempotence key) — accepted, record in DECISIONS at Task 19 — cost if wrong: provenance would not distinguish how a value was set, only what it is
Task 3: minor (deferred): pin a golden config_hash in the test; frozen=True is shallow (docstring overstates); parse_override does not strip key whitespace; json.dumps on non-JSON override leaf raises TypeError
Task 4: dispatched, BASE b6ba651
Task 4: complete (commits b6ba651..dbd7ece, review clean; line table 73 lines / 30 usable at R=23081, regeneration check deferred to Task 11 per plan)
Ruling: Task 5's HeaderMap.classify implements the rule syntax Task 4 wrote into header_map.yaml (`when: {field: {lt|eq|missing: value}}`, unconditional fallback rule allowed; rule names "rule:<i>"; "unclassified" only when no rule matches) — cost if wrong: one YAML edit
Ruling: load_anchors returns pairs only; Task 13 reads the anchors file's `trace` key itself (extend registry with a small helper if needed) — cost: one helper
Task 4: minor (deferred): config.yaml line >100 chars and crop comment wording; survey45 generator state across arc slots; untested registry error paths; load_layout/load_line_table re-validate the instrument each call; CSV header comment wrapped mid-phrase; header map ships an unconditional fallback rule (so "unclassified" is unreachable for HighSpec)
Task 5: dispatched, BASE dbd7ece
Task 5: complete (commits dbd7ece..01e3452, review clean)
Task 5: minor (deferred): unconvertible first key stops resolution (no fallthrough to next key/regex/default, no log); unmatched manifest rows dropped silently; missing folder yields [] silently; `lt` on untyped string raises TypeError; MANIFEST_FIELDS hand-duplicates FrameMeta fields; unknown top-level map keys ignored; test gaps (Undefined card, BinTableHDU skipped)
Task 6: dispatched, BASE 01e3452
Task 6: review → Needs fixes (2 Important: background frame slot; no on-chip mask for off-chip catalogue lines). FIX_BASE 3d265e2
Ruling: the simulator's drift law is anchored at the full arc: shift(t)/stretch(t) interpolate from drift endpoints[0] at the full arc's mid-time to endpoints[1] at the last frame, so the full arc's row is exactly consistent with the law (no special-casing) — resolves the brief's internal conflict — cost if wrong: one formula
Task 6: minor (deferred): shift/stretch recorded on non-drifting frames (NaN would be clearer); arc light scaled by throughput (document in SimTruth docstring); test helper cropped() ignores orientation; _fiber_fluxes returns None on unknown type; thin saturation margin (48k of 65535) — add a max<saturation assertion; blue 232 px of chip hold no catalogue line (note for wavecal tests); simulator.py 637 lines with truth I/O inside
Task 6: fix round 1/5 (2 addressed + ruling applied, 0 open; commits 3d265e2..ae1ec55)
Task 6: complete (commits 01e3452..ae1ec55, review clean after round 1)
Task 7: dispatched, BASE ae1ec55
Task 7: implementer DONE_WITH_CONCERNS at 4132fbf (square 5x5 cosmic footprint flags ~235 sky emission-line cores per simulated science frame; (1,5) leaves 4 with equal recovery)
Ruling: cosmic-ray median filter footprint defaults to along-dispersion (1, size) so emission-line cores (extended across rows) are not flagged; CalibParams.cosmic_size accepts int or (rows, cols); multi-pixel tracks along a row are a known gap until Horne profile rejection (Plan 2) / a Laplacian detector (Plan 3) — cost if wrong: some along-row cosmic rays survive into aperture extraction, where the mask would have caught them
Ruling: variance uses signal_e = gain·(raw − bias) (Poisson of the full detected signal incl. dark/background), as dispatched — cost: none (conservative)
Task 7: review → Needs fixes (2 Important: NaN/masked pixels feed the median filter unchanged; build_master averages differing exposure times silently). FIX_BASE 38845f8
Task 7: minor (deferred): median-stack variance non-robust when clip_sigma is None; MAD=0 floor drops a frame in tiny quantized stacks; footprint() rejects np.integer; no shape check across masters; all-NaN pixel variance 0.0; "2" in steps[1] weak assertion; type: ignore on shape; no-bias+background path approximates
Task 7: fix round 1/5 (2 addressed, 0 open; commits 38845f8..10cdc1e)
Task 7: complete (commits ae1ec55..10cdc1e, review clean after round 1)
Task 7: minor (deferred): median filter runs twice on clean frames; None exposure times skip the spread check; non-finite variance pixel silently never flagged
Note for Task 18 runner: group dark/background frames by exposure time before build_master (raises on mixed exposures)
Task 8: dispatched, BASE 10cdc1e
Task 8: implementer DONE_WITH_CONCERNS at 80dc470 (noise estimator departs from the notebook: statistical noise of the band median instead of the profile MAD, since the MAD of a dense 45-trace comb is the comb; classify_arc uses the mean; find_traces on a line-only arc finds nothing with a median collapse; thin calfiber margin on the dense layout)
Ruling: TraceParams gains `collapse: Literal["median","mean","percentile"] = "median"` plus `collapse_percentile: float = 90.0`; the arc_full trace source uses the percentile collapse so line-only frames yield a trace profile; the runner passes it when trace_source returns "arc_full" — cost if wrong: one parameter
Ruling: the collapsed-profile noise estimator = statistical noise of the band median (not the profile MAD) is accepted; the reviewer must confirm lab3 (3 traces) and survey45 (45) simulated flats both give exactly the layout's count — cost if wrong: Task 19's lab regression will show it
Task 8: review → Needs fixes (3 Important: classify_arc absolute threshold fragile; match_layout can alias by one pitch on the periodic survey layout; comb_shift kernel width coupled to step). FIX_BASE 30e53ae
Task 8: minor (deferred): _clipped_fit may return coeffs from before the last clip; collapse_percentile range-check at runtime not in the model; redundant NaN test; np.std clip scale; trailing columns never centroided; blanket RuntimeWarning ignore in _collapse; weak `< 45` assertion; noise==0 edge in find_peaks height
Task 8: fix round 1/5 (3 addressed, 0 open; commits 30e53ae..f23a5ec)
Task 8: complete (commits 10cdc1e..f23a5ec, review clean after round 1)
Task 8: minor (deferred): comb_shift width test does not discriminate; brightest level in classify_arc includes unmatched traces; lit_fraction=0.1 makes a full arc with a <10% fiber come back unknown (document when surfaced in config)
Task 9: dispatched, BASE f23a5ec
Task 9: implementer DONE_WITH_CONCERNS at 0e1ed9a (observations: calibrate_frame variance omits the background master's noise — pull sigma 1.29 vs 1.28 predicted; trace._aperture_sum not unified with extract (circular import); relative_throughput uses the other aperture)
Task 7: minor (deferred, found in Task 9): calibrate_frame should propagate MasterCalib.variance of the subtracted masters into the frame variance (currently ~1% error underestimate on the simulator)
Task 9: complete (commits f23a5ec..0e1ed9a, review clean; ⚠️ calibrate_frame marks non-finite in mask — yes, Task 7 masks non-finite values)
Task 9: minor (deferred): NaN pixels reduce coverage but never flag the column when coverage stays above min_coverage; n_iter semantics off by one; vacuous mask assertion in test; ~1e-16 subpixel edge weights; robust_sigma should move to a shared helper to let trace._aperture_sum reuse _aperture_grid
Task 10: dispatched, BASE 0e1ed9a
Task 10: implementer DONE_WITH_CONCERNS at 4368496 (use_y not barred for saturated lines; _numeric_fwhm NaN guard departs from notebook; narrow deblend component at sigma bound rejected only by on_trace; 16 s module fixture)
Task 10: review → Needs fixes (2 Important: `complete` in use_x computed over the full strip voids lines on traces within halfwidth of the detector edge; deblender can fabricate a component pinned at the sigma lower bound). FIX_BASE 4368496
Ruling: use_x follows the brief (fit_ok|isolated, ~near_edge, on_trace, ~saturated), `complete` is dropped from use_x and computed on the measuring rows only; the deblender rejects a trial whose new component's sigma sits at sigma_bounds[0] or below 0.5·sigma_guess — cost if wrong: a few lines excluded/included at the margins
Task 10: minor (deferred): x_prof recomputed per line; GAUSSIAN_FWHM duplicated between lines.py and simulator; find_lines ~230 lines; _fit_group catches ValueError broadly; weakest<add_sigma test applies to all components; pass-1 peak/component pairing by index; calibrated-full-arc fixture duplicated in three test modules (make a session conftest fixture)
Task 10: fix round 1/5 (2 addressed, 0 open; commits 4368496..f18f7c5)
Task 10: complete (commits 0e1ed9a..f18f7c5, review clean after round 1)
Task 10: minor (deferred): "newly added component" identified by nearest-to-residual-peak not by identity; CHANGES.md lacks the two departures; initial _fit_group can still return a sigma-bound component (on_trace's job); use_y not barred for saturated lines
Task 11: dispatched, BASE f18f7c5
Task 11: implementer DONE at a89ab0b; regeneration check changed the line table by one threshold-straddling line (30→31 usable; 514.6055 Th I at SNR 20.44 with the notebook's truncated np.convolve kernel vs 19.996 with gaussian_filter1d). Ruling: the notebook-faithful kernel wins (31 is the notebook's own count); regenerated table committed — cost if wrong: one marginal line in the catalogue
Note for Task 12: fit_model starts the distortion term at zero every call; warm-start only x_c, B, alpha from a previous reference
Task 11: complete (commits f18f7c5..a89ab0b, review clean)
Task 11: minor (deferred): grating_x raises inside least_squares for wild alpha steps (guard or alpha_bound default); fit_model discards a start model's poly when poly_deg is below the lowest power; Protocol lacks from_dict; atlas noise==0 division; CRPIX1 assumption undocumented; duplicated unsupported-kind branches
Task 12: dispatched, BASE a89ab0b
Task 12: implementer DONE at 32c1365 (45/45 solved, rms 0.005-0.015 px, max |dλ|/disp 0.06; solve_full_arc(..., ref_fiber=0, geom=None) geom required for a ReferenceSolution start; conftest gained session fixtures sim_full_arc/sim_arc_lines/sim_arc_fluxes/truth_traceset). Note for Task 17/18: the geometry must be rebuilt from the instrument config when reading a reference
Task 12: complete (commits a89ab0b..32c1365, review clean)
Task 12: minor (deferred): warm-start test does not assert seed provenance; `<= {3,5}` passes on empty; vacuous hashlib assert; silent zero-shift fallback on reference/traceset mismatch; empty fiber_agreement types n_fibers float; give-up table dx_px semantics; single-model start skips the n_pix check; zero inter-fiber shift in the simulator leaves the shift-seeding paths weakly tested; test_lines imports from tests.conftest
Task 13: dispatched, BASE 32c1365
Ruling: Task 13 (identify) writes the reference FITS product defined in Task 17, a forward dependency the pre-flight scan missed; every dataclass Task 17 serialises already exists (Task 2), so Task 17 is executed before Task 13 (order 17, 13, 14, 15, 16, 18, 19) — cost if wrong: none (Task 17 has no dependency on 13–16 beyond the dataclasses)
Task 13: dispatch withdrawn (not started); Task 17: dispatched, BASE 32c1365
Task 17: implementer DONE at b2f6b1b (read_spectra → (es, prov, drift|None); read_combined → (cs, meta, prov); write_spectra(..., drift=, drift_scatter_px=, drift_usable=); read_reference re-hashes; non-ASCII targets unhandled)
Task 17: review → Needs fixes (2 Important: write_spectra silently overrides Provenance statuses/reference/drift arcs with the ExposureSpectra's; drift scatter/usable carried in DataFrame.attrs). FIX_BASE b2f6b1b
Ruling: write_spectra takes the drift as a `DriftModel` (its per_fiber table plus scatter_px/usable) with an optional `weights: dict[int,float] | None`; read_spectra returns `(es, prov, DriftModel | None)`; provenance fields that both es and prov carry must agree or write_spectra raises ValueError — cost if wrong: one signature change for Task 18
Task 17: minor (deferred): presence card wins over HAS_* rows; _check_geometry hasattr-skips; DRIFT written independent of drift_arcs; untested error paths; DRIFT header not in the raw-keyword test; ""↔None conflation; _tuplify lists; _set relies on astropy card widths; NaN MPSHIFT; products.py 914 lines
Task 17: fix round 1/5 (2 addressed, 0 open; commits b2f6b1b..c984fdb)
Task 17: complete (commits 32c1365..c984fdb, review clean after round 1)
Task 17: minor (deferred): Provenance flux/telluric defaults count as "set" (asymmetric unset rule) — Plan 2's flux-cal runner must pass the exposure's status explicitly; WEIGHT column written but not surfaced by read_spectra
Note for Task 18: write_spectra(path, es, prov, drift: DriftModel|None, weights: dict|None); read_spectra → (es, prov, DriftModel|None); read_combined → (cs, meta, prov); read_reference(path, geometry_lookup)
Task 13: dispatched, BASE c984fdb
Task 13: implementer DONE at 8ee9cab (find_lines ~30 s per 45-trace chain dominates; suite ~3m20s; run_identify/run_solve take the layout preset name, _arc_lines takes the table)
Task 13: review → Approved with 2 Important (provenance lists a background master even when subtract_background is off; _latest_master sorts by name so notime/other-prefix masters win) → fix loop. FIX_BASE 8ee9cab
Task 13: minor (deferred): _arc_lines hard-wires out_dir/calib globs (Task 18 runner needs a masters/traceset seam); _interp empty → opaque ValueError; synthetic template keeps pedestal; fiber_agreement.csv fixed name overwritten; read_anchors_file duplicates registry pair parsing; assert as control flow; anchors=[] sentinel; trace index is detection order (layout_index would be more robust); identify.py 555 lines; slow marker for survey-night tests
Task 13: fix round 1/5 (2 Important + 1 Minor seam addressed, 0 open; commits 8ee9cab..b1d8da8)
Task 13: complete (commits c984fdb..b1d8da8, review clean after round 1)
Task 13: minor (deferred): DECISIONS.md still says masters are "latest by name"; background "saturation mask only" step over-states when detector.saturation is None; noconf masters never match a named configuration
Task 14: dispatched, BASE b1d8da8
Task 14: implementer DONE at 0ff2d46 (sim cal arc: shift error ≤0.0011 px, transfer ≤0.136 disp; comb_shift argument order swapped to get the drift sense; warn_shift_px left to the runner)
Task 14: complete (commits b1d8da8..0ff2d46, review clean)
Ruling: comb_shift(x_ref, x_obs) argument order in the drift fit (returns the drift-sense shift; the brief's literal order gives the negative) — accepted — cost: none
Task 14: minor (deferred): apply_drift docstring lacks the weights-as-given contract; None configurations pair with each other; drift tests import private _arc_lines; rank-deficient polyfit with identical y_centers; transfer test tolerance absorbs arc→science drift
Task 15: dispatched, BASE 0ff2d46
Task 15: implementer DONE at faf2cce (linear resampling leaves up to ~7% of a 3-px line peak between sub-pixel-offset grids — note for Plan 2/3; truth_traceset fixture now sets telescope)
Task 15: review → Needs fixes (1 Important: resampled sky error = interpolated σ, up to 41% high at half-pixel offsets; exact = sqrt((1−t)²σ0² + t²σ1²)). FIX_BASE faf2cce
Task 15: minor (deferred): NaN/non-positive throughput unguarded (science → NaN skysub unmasked; sky → silent ratio 1.0); RESAMPLERS duplicates the Literal; duplicate fiber in by_fiber silently last-wins; conftest isinstance-float NaN test; night fixture duplicates calibration chain (make a shared science-frame fixture); aliasing in a test
Task 15: fix round 1/5 (1 Important + 1 Minor addressed, 0 open; commits faf2cce..f2d136a)
Task 15: complete (commits 0ff2d46..f2d136a, review clean after round 1)
Task 15: minor (deferred): single-pixel sky spectrum hits 0/0 under -W error; None throughput entry raises TypeError instead of "bad throughput"
Task 16: dispatched, BASE f2d136a
Task 16: implementer DONE at 037173d (clipping guarded so a pixel is never emptied; resample_sky t≈1e-13 edge noted; source column empty)
Task 16: review → Needs fixes (1 Important: clipping deviation divides by error==0 → RuntimeWarning under -W error). FIX_BASE 037173d
Task 16: minor (deferred): fully-masked contributor listed with weight 0 instead of excluded; _reject keyed by fiber; telluric imports the helper from fluxcal; misnamed test; _coverage median over non-finite steps
Task 16: fix round 1/5 (1 Important + 2 requested items addressed, 0 open; commits 037173d..93ecbed)
Task 16: complete (commits f2d136a..93ecbed, review clean after round 1)
Ruling: pytest config gets `filterwarnings = ["error"]` so the pristine-output constraint is enforced in CI, not only in manual `-W error` runs — folded into Task 18 — cost: a warning anywhere fails CI (intended)
Task 18: dispatched, BASE 93ecbed
Task 18: implementer DONE_WITH_CONCERNS at 9786848 (suite 222 passed/1 skipped in 6m29s under filterwarnings=error; reduced_again covers any key change; identify/solve register in the out folder's CalibLibrary; _arc_lines accepts a FrameRef; byte-truncated FITS only warns in astropy)
Plan 2 note: solutions/drift/sky are keyed by detection-order trace index; a missing trace misaligns reference fibers — key by layout_index in Plan 2
Task 18: review → Needs fixes (2 Important: drift models stored without the reference they were fitted against, so a multi-reference night under nearest/interpolate can apply B-relative drift to A's wavelengths; the interpolate blend path (_blend_drifts) is untested). FIX_BASE 9786848
Task 18: minor (deferred): WEIGHT column holds summed normalised weights (blend split only implicit); pinned --reference not in library → per-frame KeyErrors instead of one group failure; _timestamp duplicated in callib and arcs; sorted(groups.items(), key=str); `if flat_frame` truthiness; registry_geometry returns Any; --layout not in the idempotence key; solve --library takes masters from out/calib; masters rebuilt every run; loose drift_shift atol in e2e test; runner.py 1075 lines; show renderers could move out of commands.py
Task 18: fix round 1/5 (2 Important + 1 Minor addressed, 0 open; commits 9786848..e6b6a51)
Task 18: complete (commits 93ecbed..e6b6a51, review clean after round 1)
Task 18: minor (deferred): _arc_weights drops an absent checksum silently; e2e test requests an unused fixture
Task 19: dispatched, BASE e6b6a51
Task 19: implementer DONE_WITH_CONCERNS at f6ea1ec (regression vs the notebook: traces 17.51/188.96/360.62 px, solution rms 0.016/0.033/0.012 px, lines 23/25/20, coverage to 0.03 nm, agreement median 0.023 px, self-drift 0.0; suite 228 passed/5 skipped without data, 233/0 with)
Ruling (for the final fix wave): when a layout is given, the runner extracts only traces matched to the layout (role not None) and logs the `extra` detections (real flat shows two faint ambient-light features at 5σ) — cost if wrong: an unlisted real fiber would be dropped, but the layout is the authority by design
Ruling (for the final fix wave): header-map field entries gain an optional `missing_values` list; HighSpec `detector_temperature` lists -999 (Andor "unstabilised") so it resolves to None instead of a real temperature — cost: one YAML key
Task 19: complete (commits e6b6a51..f6ea1ec, review clean)
Task 19: minor (deferred): README "every command takes --instrument" untrue of show; --reference and --anchors-file flags not listed; manifest field list omits classification_rule; slow-marker description; brittle exact-zero stretch assertion
All 19 tasks complete. Final whole-branch review: MERGE_BASE dd1f8d6, HEAD f6ea1ec
Final review (dd1f8d6..f6ea1ec): With fixes — 1 Critical (per-fiber keys follow detection-order trace index; a missing/extra trace between nights mis-assigns solutions silently), 6 Important (runner extracts unmatched traces; TEMP=-999 as real temperature; mid_time is start time; survey45 defaults overlap fibers; missing/empty night folder exits 0; reference/frame n_pix mismatch caught late), minors; ledger triage given. FIX_BASE f6ea1ec
Ruling (reversal): re-key fibers by layout_index now, not in Plan 2 — the final reviewer showed the wrong-science path and the ingredients exist — cost if wrong: a wide refactor late in the plan, mitigated by a new regression-style test
Ruling: CLAUDE.md edits (lint paths, phantom `reduce`/`focus` commands) need the user's approval — left out of the fix wave, surfaced to the user
Ruling: barycentric corrections (spec requires; no task scheduled them; no pointing in the header vocabulary) → Plan 2 — cost: CombinedSpectrum.barycentric stays {} in Plan 1 products
Final fix wave: dispatched (one subagent, complete findings list)
Final fix wave: DONE at e5abbbb (240 passed/5 skipped; 245/0 with lab data; regression unchanged within tolerances; apply_drift now returns (wavelengths, skipped); solve_full_arc fibers/ref_fiber are trace ids; mid-time derived so real-arc reference stamps shift by half the exposure)
