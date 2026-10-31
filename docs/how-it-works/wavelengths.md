# Wavelength calibration in three cadences

Arc light reaches the fibers in two very different ways, and the calibration
follows that. Feeding the twenty telescope fibers from a lamp is expensive and
happens every few months; the five calibration fibers can be lit from the
calibration unit in the instrument room at any time, and are, one to three
times a night and at every instrument switch. Science exposures carry no arc
light at all.

## 1. The full arc: a new reference solution

Every few months every trace sees the lamp. Each fiber gets a **fresh,
independent** solution: Th-Ar lines are found on the trace, matched to the
catalogue through a starting solution, and a grating-equation model is fitted
and refined in three passes of identify → fit → clip until the wavelength map
stops moving. The starting solution is the previous reference for that fiber
(the first time, two anchor pairs read off the panels by hand; see
[The first reference solution](../running/first-reference.md)); it is only a
seed and carries the line identities forward, it does not constrain the fit.

The model is the notebook's: the grating equation with the incidence angle held at
its design value (it is degenerate with the plate scale at 0.9993 correlation)
and a low-order distortion polynomial. The anchors fix the two parameters of the pure grating model, the pixel of
the Littrow wavelength and the plate scale; the full fit adds the low-order
distortion terms. The curvature of λ(x) comes from the physics rather than
from a polynomial. On the laboratory
frame the three fibers solve with 23, 25 and 20 lines at 0.016, 0.033 and
0.012 px rms, within the notebook's numbers; because the fibers are solved
independently, the agreement between them on common lines (median 0.02 px) is
an external accuracy check, not a fit statistic.

The result is a **reference solution**: per fiber the model, the identified
lines with residuals, the rms, the pixel and wavelength range the lines span,
and the trace centre; plus an identifier, a hash, and a validity that starts
at the arc's mid-time and ends at the next full arc of the same grating and
layout. References are never edited; a new full arc makes a new one.

## 2. The calibration-fiber arc: a drift model

Between full arcs only the five calibration fibers speak. Each of their arcs
is reduced to lines; for each arc fiber the measured line positions are
matched to the reference solution's identified lines and a **linear drift in
pixel space** is fitted:

```text
x_observed − x_reference = shift + stretch · (x_reference − n/2)
```

with the stretch set to zero when it is not significant. A low-order
polynomial in trace position smooths shift and stretch across the detector;
its scatter about the five fibers is the drift model's quality number. On the
simulated night the planted shift is recovered to well within the 0.05 px the
tests demand.

## 3. The science frame: transfer by time

A science frame takes its wavelengths from the reference valid at its
mid-time, corrected by the drift model of the calibration-fiber arcs **paired
by time**: the latest before it (`nearest_before`, the quick-look default),
the nearest either way, or the arcs before and after weighted by time
(`interpolate`, for the manual reduction once the night is complete). Only
arcs of the same grating, fitted against the same reference, and within
`max_age_hours` qualify. For each fiber the drift evaluated at its trace
position is undone and the reference model gives the wavelength:

```text
x_reference = (x − shift + stretch · n/2) / (1 + stretch)
λ(x) = reference_model(x_reference)
```

The frame's **status** records what happened: `transferred`, or `reference`
when no arc paired, or `nominal` (geometry only) when no reference exists.
The product also names the reference identifier and hash, the drift arcs'
checksums and weights, and the fitted drift per fiber, so a spectrum reduced a
month into a run is traceable to the full arc it rests on.

## Why not solve every night

The notebook solved each frame from scratch with a by-eye anchor. That is the
right thing for a laboratory frame and the wrong thing for an unattended
service: the hard, attended part now happens once per grating and layout (`identify`),
the full arcs are solved automatically from the previous reference, and the
nightly job is a two-parameter fit per calibration fiber. The design also
matches the hardware: science fibers cannot be lit nightly, so a nightly full
solve was never available.

## Fiber identity

Everything per fiber is keyed by the fiber layout's trace id, not by the order
in which traces were detected. A fiber that is missing on a night is absent
from its products; a trace whose centre has moved by more than the layout
tolerance from the reference solution's is given no wavelength and listed in
the QC. Without this, a missing or extra trace between the reference night and
a science night would hand each fiber its neighbour's solution.

## Air, vacuum and the barycentre

Wavelengths are **air wavelengths in nanometres** throughout, matching the
Murphy et al. atlas. Barycentric corrections are to be computed per telescope
and stored, never applied to the grid; they are not yet computed in this
release because the header vocabulary carries no pointing.
