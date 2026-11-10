"""Stages 4 and 5 of the quick-look recipe: a science frame calibrated, extracted, given
wavelengths by the time policy, sky-subtracted and written with its provenance and QC,
then combined per target. The whole recipe is described in :mod:`mastpipe.pipeline.runner`.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.__version__ import __version__
from mastpipe.core.model import PER_FIBER_COLUMNS, DriftModel, ExposureSpectra, FrameMeta, FrameRef, ReferenceSolution
from mastpipe.core.times import to_utc, utc_iso
from mastpipe.io import products
from mastpipe.pipeline.arcs import pair_arcs, status_for
from mastpipe.pipeline.group import (
    ALREADY_REDUCED,
    _calibrated,
    _check_binning,
    _check_readout,
    _configuration_origin,
    _configuration_steps,
    _content_steps,
    _describe,
    _Group,
    _median_snr,
    _reference_for,
    _skip,
    _write_qc,
    log,
)
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.combine import combine_spectra, group_by_target
from mastpipe.steps.drift import apply_drift, drift_unusable_reason
from mastpipe.steps.extract import extract_aperture
from mastpipe.steps.fluxcal import apply_response, flux_status
from mastpipe.steps.sky import subtract_sky_all
from mastpipe.steps.telluric import correct_telluric, telluric_status
from mastpipe.steps.wavemodel import model_from_config


def _arcs_for(
    arcs: list[tuple[FrameMeta, DriftModel, str]], reference: ReferenceSolution
) -> list[tuple[FrameMeta, DriftModel]]:
    """The drift models fitted against ``reference``: a shift measured against another
    reference (a full arc solved mid-night) must never be applied to this one's map."""
    return [(meta, drift) for meta, drift, identifier in arcs if identifier == reference.identifier]


def _arc_weights(pairs: list[tuple[DriftModel, float]], checksums: list[str]) -> list[float]:
    """The pairing weight of every checksum of ``checksums`` (the arcs' union), for QC."""
    weights = {c: 0.0 for c in checksums}
    for model, weight in pairs:
        for checksum in model.arc_checksums:
            weights[checksum] = weights.get(checksum, 0.0) + float(weight)
    return [weights[c] for c in checksums]


def _blend_drifts(
    pairs: list[tuple[DriftModel, float]],
) -> tuple[DriftModel, dict[int, float] | None]:
    """One :class:`~mastpipe.core.model.DriftModel` for the product: the model itself for a
    single arc (weights ``None``: every fiber 1.0); for several, the weighted blend of their
    polynomials (right-aligned when the degrees differ) and per-fiber values (``n_lines`` the
    smallest, ``usable`` when all are), the sorted union of their checksums, the weighted
    mid-time, the largest scatter, the smallest ``n_usable`` and ``usable`` when every model
    is. The returned weights feed the product's ``WEIGHT`` column: per fiber, the summed
    normalised weight of the models that hold the fiber (1.0 when all do), which says how
    much of the blend that row represents; the per-arc weights themselves are recorded in
    the science QC JSON (``drift_weights``, aligned with ``drift_arcs``)."""
    if len(pairs) == 1:
        return pairs[0][0], None
    total = sum(weight for _, weight in pairs)
    pairs = [(model, weight / total) for model, weight in pairs]
    n_shift = max(len(model.shift_coeffs) for model, _ in pairs)
    n_stretch = max(len(model.stretch_coeffs) for model, _ in pairs)
    shift = np.zeros(n_shift)
    stretch = np.zeros(n_stretch)
    for model, weight in pairs:
        shift[n_shift - len(model.shift_coeffs) :] += weight * np.asarray(model.shift_coeffs)
        stretch[n_stretch - len(model.stretch_coeffs) :] += weight * np.asarray(model.stretch_coeffs)
    tables = [(model.per_fiber.set_index("fiber"), weight) for model, weight in pairs]
    fibers = sorted({int(f) for table, _ in tables for f in table.index})
    rows: list[dict[str, Any]] = []
    fiber_weights: dict[int, float] = {}
    for fiber in fibers:
        present = [(table.loc[fiber], weight) for table, weight in tables if fiber in table.index]
        weight_sum = sum(weight for _, weight in present)
        fiber_weights[fiber] = float(weight_sum)
        rows.append(
            {
                "fiber": fiber,
                "y_center": float(present[0][0]["y_center"]),
                "shift_px": sum(w * float(r["shift_px"]) for r, w in present) / weight_sum,
                "stretch": sum(w * float(r["stretch"]) for r, w in present) / weight_sum,
                "n_lines": int(min(int(r["n_lines"]) for r, _ in present)),
                "rms_px": sum(w * float(r["rms_px"]) for r, w in present) / weight_sum,
                "usable": bool(all(bool(r["usable"]) for r, _ in present)),
            }
        )
    origin = to_utc(pairs[0][0].arc_mid_time)
    offset = sum(
        ((to_utc(model.arc_mid_time) - origin).total_seconds() * weight for model, weight in pairs),
        0.0,
    )
    mid_time = utc_iso(origin + pd.Timedelta(seconds=offset))
    blended = DriftModel(
        arc_checksums=tuple(sorted({c for model, _ in pairs for c in model.arc_checksums})),
        arc_mid_time=mid_time,
        per_fiber=pd.DataFrame(rows, columns=list(PER_FIBER_COLUMNS)),
        shift_coeffs=shift,
        stretch_coeffs=stretch,
        scatter_px=max(float(model.scatter_px) for model, _ in pairs),
        usable=all(model.usable for model, _ in pairs),
        n_pix=pairs[0][0].n_pix,
        n_usable=min(int(model.n_usable or 0) for model, _ in pairs),
    )
    return blended, fiber_weights


def _wavelengths(
    group: _Group,
    reference: ReferenceSolution | None,
    pairs: list[tuple[DriftModel, float]],
    n_pix: int,
) -> tuple[str, dict[int, np.ndarray], DriftModel | None, dict[int, float] | None, dict[int, str]]:
    """``(status, wavelength by fiber, drift model, weights, skipped)`` of a science frame.

    ``skipped`` maps the fibers that keep no wavelength to the reason (a trace centre
    beyond ``trace.layout_tolerance`` from the reference solution's, see
    :func:`~mastpipe.steps.drift.apply_drift`). ``ValueError`` when the reference was
    solved on another pixel grid than the frame's ``n_pix`` columns.
    """
    if reference is None:
        model = model_from_config(group.effective.geometry.model, group.geom, n_pix)
        lam = np.asarray(model.wavelength(np.arange(n_pix, dtype=float)), dtype=float)
        return "nominal", {trace.index: lam.copy() for trace in group.traces()}, None, None, {}
    grids = sorted({int(solution.n_pix) for solution in reference.solutions.values()})
    if grids and grids != [int(n_pix)]:
        raise ValueError(
            f"reference {reference.identifier} is solved on {grids} pixels, the frame has "
            f"{n_pix} columns: the reference belongs to another detector crop"
        )
    wavelengths, skipped = apply_drift(
        reference,
        pairs,
        group.traces(),
        group.geom,
        tolerance=group.params.trace.layout_tolerance,
        n_pix=n_pix,
    )
    if not pairs:
        return status_for(pairs), wavelengths, None, None, skipped
    drift, weights = _blend_drifts(pairs)
    return status_for(pairs), wavelengths, drift, weights, skipped


def _idempotence_key(
    checksum: str | None,
    config_hash: str,
    reference_hash: str | None,
    drift_arcs: Iterable[str],
    layout: str | None,
    layout_hash: str | None,
    masters_hash: str | None,
    timezone: str | None,
) -> tuple[Any, ...]:
    """A science product's idempotence key (described in :mod:`mastpipe.pipeline.runner`),
    the drift-arc checksums sorted: built alike for the run (:func:`_science_frame`) and
    for the product on disk (:func:`_stored_key`), so the two compare equal exactly when
    every part does."""
    return (
        checksum,
        config_hash,
        reference_hash,
        tuple(sorted(drift_arcs)),
        layout,
        layout_hash,
        masters_hash,
        timezone,
    )


def _stored_key(path: Path) -> tuple[tuple[Any, ...], str] | None:
    """The idempotence key and wavelength status of the product at ``path``; ``None`` when
    there is none (an unreadable product counts as none and is rewritten)."""
    if not path.is_file():
        return None
    try:
        prov = Provenance.from_header(fits.getheader(path))
    except (OSError, KeyError, ValueError) as error:
        log.warning("%s cannot be read as a product (%s); it will be rewritten", path.name, error)
        return None
    key = _idempotence_key(
        prov.inputs[0] if prov.inputs else None,
        prov.config_hash,
        prov.reference_hash,
        prov.drift_arcs,
        prov.layout,
        prov.layout_hash,
        prov.masters_hash,
        prov.timezone,
    )
    return key, prov.wavelength_status


def _science_frame(group: _Group, ref: FrameRef) -> None:
    """Stage 4 for one science frame: wavelengths by the time policy, sky, products, QC.
    A frame not binned or read out as the group
    (:func:`~mastpipe.pipeline.group._check_binning`,
    :func:`~mastpipe.pipeline.group._check_readout`) fails rather than being calibrated
    with masters and noise of another binning or amplifier."""
    _check_binning(group, ref)
    _check_readout(group, ref)
    meta, params = ref.meta, group.params
    reference = _reference_for(group, meta.mid_time)
    arcs = [] if reference is None else _arcs_for(group.arcs, reference)
    pairs = [] if reference is None else pair_arcs(meta, arcs, reference, group.drift_params)
    drift_arcs = sorted({c for model, _ in pairs for c in model.arc_checksums})
    reference_hash = None if reference is None else (reference.hash or None)
    key = _idempotence_key(
        meta.checksum,
        group.cfg_hash,
        reference_hash,
        drift_arcs,
        group.preset,
        group.layout_hash,
        group.masters_hash,
        meta.timezone,
    )
    path = group.spec_dir / products.product_name(meta, "spec")
    stored = _stored_key(path)
    if stored is not None:
        if stored[0] == key and not group.options.force:
            _skip(group, meta.source, f"{ALREADY_REDUCED} as {path.name}")
            return
        why = "forced" if stored[0] == key else f"its key changed (was {stored[1]})"
        log.info("reducing %s again: %s", meta.source, why)

    frame, steps = _calibrated(group, ref, cosmics=True)
    steps = [*_content_steps(meta), *_configuration_steps(group, meta), *steps]
    traceset = group.traces()
    spectra = extract_aperture(frame, traceset, params.extract)
    steps.append(
        f"extract {params.extract.aperture} aperture of halfwidth {params.extract.halfwidth:g} ({params.extract.method})"
    )
    status, wavelengths, drift, weights, skipped = _wavelengths(group, reference, pairs, frame.shape[1])
    for spectrum in spectra:
        if spectrum.fiber in wavelengths:
            spectrum.wavelength = wavelengths[spectrum.fiber]
    withheld = [
        {"fiber": int(fiber), "reason": reason}
        for fiber, reason in sorted(skipped.items())
        if any(spectrum.fiber == fiber for spectrum in spectra)
    ]
    for item in withheld:
        log.warning("%s: fiber %d is given no wavelength: %s", meta.source, item["fiber"], item["reason"])
    if reference is None:
        steps.append("wavelength nominal from the configuration geometry")
    elif drift is None:
        steps.append(f"wavelength from reference {reference.identifier} (no drift arc paired)")
    else:
        steps.append(
            f"wavelength transferred from reference {reference.identifier} through "
            f"{len(pairs)} drift arc(s) ({group.drift_params.time_policy})"
        )
    if withheld:
        steps.append(f"wavelength withheld from {len(withheld)} fiber(s) whose trace is off the reference's centre")
    spectra, sky_reasons = subtract_sky_all(spectra, traceset, group.throughput, params.sky)
    n_sky = sum(1 for spectrum in spectra if "skysub" in spectrum.versions)
    steps.append(f"skysub: {n_sky} science fibers subtracted, {len(sky_reasons)} left raw")
    spectra = [correct_telluric(apply_response(spectrum, None), None) for spectrum in spectra]
    es = ExposureSpectra(
        meta=meta,
        spectra=spectra,
        wavelength_status=status,
        flux_status=flux_status(spectra),
        telluric_status=telluric_status(spectra),
        reference_id=None if reference is None else reference.identifier,
        reference_hash=reference_hash,
        drift_arcs=tuple(drift_arcs),
    )
    prov = Provenance(
        pipeline_version=__version__,
        config_hash=group.cfg_hash,
        inputs=[meta.checksum],
        calibrations=group.calibration_names(),
        wavelength_status=es.wavelength_status,
        flux_status=es.flux_status,
        telluric_status=es.telluric_status,
        reference_id=es.reference_id,
        reference_hash=es.reference_hash,
        drift_arcs=list(drift_arcs),
        steps=steps,
        layout=group.preset,
        timezone=meta.timezone,
        layout_hash=group.layout_hash,
        masters_hash=group.masters_hash,
    )
    reason = None if drift is None else drift_unusable_reason(drift, group.drift_params)
    products.write_spectra(path, es, prov, drift=drift, weights=weights, drift_reason=reason)
    group.summary.products.append(path)
    if stored is not None:
        group.summary.reduced_again.append(meta.source)
    traces = {trace.index: trace for trace in traceset}
    _write_qc(
        group,
        path.with_suffix(".json").name,
        "science",
        {
            "source": meta.source,
            "classification_rule": meta.classification_rule,
            "configuration_origin": _configuration_origin(group, meta),
            "mid_time": meta.mid_time,
            "targets": dict(meta.targets),
            "wavelength_status": status,
            "reference_id": es.reference_id,
            "drift_arcs": list(drift_arcs),
            "drift_weights": _arc_weights(pairs, drift_arcs),
            "drift_scatter_px": None if drift is None else float(drift.scatter_px),
            "wavelength_skipped": withheld,
            "fibers": [s.fiber for s in spectra],
            "roles": [s.role for s in spectra],
            "telescopes": [s.telescope for s in spectra],
            "median_snr": [_median_snr(s) for s in spectra],
            "masked_fraction": [float(np.mean(s.mask)) for s in spectra],
            "drift_shift_px": [None if drift is None else float(drift.shift_at(traces[s.fiber].y_center)) for s in spectra],
            "drift_stretch": [None if drift is None else float(drift.stretch_at(traces[s.fiber].y_center)) for s in spectra],
            "sky": sky_reasons,
        },
    )
    log.info(
        "science %s: %d spectra, wavelength %s, %d sky-subtracted -> %s",
        meta.source,
        len(spectra),
        status,
        n_sky,
        path.name,
    )
    _combine_targets(group, es, prov)


def _combine_targets(group: _Group, es: ExposureSpectra, prov: Provenance) -> None:
    """Stage 5: one combined product per target of the exposure (fibers without a target
    are not combined)."""
    meta = es.meta
    for target, members in group_by_target(es.spectra).items():
        if any(spectrum.target is None for spectrum in members):
            log.info("%s: %d science fibers without a target are not combined", meta.source, len(members))
            continue
        try:
            cs = combine_spectra(es.spectra, target, group.params.combine, version=None)
            path = group.spec_dir / products.combined_name(meta, target)
            steps = [*prov.steps, f"combine {len(cs.contributors)} science fibers of {target}"]
            products.write_combined(path, cs, meta, replace(prov, steps=steps))
            group.summary.products.append(path)
            log.info(
                "combined %s: %d fibers, %d excluded -> %s",
                target,
                len(cs.contributors),
                len(cs.excluded),
                path.name,
            )
        except Exception as error:
            if group.options.on_error == "raise":
                raise
            group.summary.failed.append((meta.source, f"combine {target}: {_describe(error)}"))
            log.exception("failed combining %s of %s", target, meta.source)
