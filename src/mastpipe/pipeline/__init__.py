"""Orchestration layer: what ties steps and I/O into commands.

:mod:`mastpipe.pipeline.provenance` is the record every product header carries and
:mod:`mastpipe.pipeline.identify` the attended bootstrap of a configuration's first
reference solution (``mastpipe identify``) with the seeded re-solve behind ``solve``. The
calibration library, arc pairing and the quick-look runner join this package in later
tasks.
"""
