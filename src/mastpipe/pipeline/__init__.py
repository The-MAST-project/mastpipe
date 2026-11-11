"""Orchestration layer: what ties steps and I/O into commands.

:mod:`mastpipe.pipeline.provenance` is the record every product header carries,
:mod:`mastpipe.pipeline.identify` the attended bootstrap of a configuration's first
reference solution (``mastpipe identify``) with the seeded re-solve behind ``solve``,
:mod:`mastpipe.pipeline.callib` the calibration library of masters and time-keyed
references, :mod:`mastpipe.pipeline.arcs` the time pairing of exposures with drift arcs,
:mod:`mastpipe.pipeline.runner` the quick-look recipe (``mastpipe quicklook``), whose
stages are :mod:`mastpipe.pipeline.masters`, :mod:`mastpipe.pipeline.traces`,
:mod:`mastpipe.pipeline.arc_stage` and :mod:`mastpipe.pipeline.science` over the shared
group state and helpers of :mod:`mastpipe.pipeline.group`, and
:mod:`mastpipe.pipeline.commands` the ``quicklook``, ``solve`` and ``show`` commands.
"""
