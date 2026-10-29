# Changes

## 0.0.1 (unreleased)

- Package skeleton: `pyproject.toml` with the `src` layout, `mastpipe` console
  script (`--version`, help), CI matrix (Linux, Windows, macOS × Python 3.11,
  3.12) and the privacy scan for tracked files.
- Core data model (`mastpipe.core.model`): the dataclasses exchanged between
  layers (`FrameMeta`, `Frame`, `MasterCalib`, `CalibrationSet`, `Trace`,
  `TraceSet`, `Spectrum` with kept `versions`, `ExposureSpectra`,
  `WavelengthSolution`, `ReferenceSolution`, `DriftModel`, `CombinedSpectrum`,
  `QCReport`) and the vocabulary constants they validate against.
- Configuration layer (`mastpipe.config`): frozen Pydantic models that forbid unknown
  keys (`DetectorConfig`, `GeometryConfig`, the per-step `*Params` under `StepParams`,
  `ConfigurationBlock`, `InstrumentConfig`, `EffectiveConfig`) and the loader
  (`load_yaml`, `deep_merge`, `build_effective` layering instrument → configuration
  block → overrides, `parse_override` for `key.path=value`, `config_hash`).
