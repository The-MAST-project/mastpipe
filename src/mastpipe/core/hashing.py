"""The length of the pipeline's short content hashes, shared by every module that makes
one. This module imports nothing, so any layer may import it without a cycle."""

HASH_LENGTH = 12
"""Hex characters kept of a SHA-256 digest in a short hash: the configuration and
calibration hashes (:func:`~mastpipe.config.loader.config_hash`,
:func:`~mastpipe.config.loader.calib_hash`), the layout hash
(:func:`~mastpipe.instruments.registry.layout_hash`), the hash of the masters' checksums
(:func:`~mastpipe.pipeline.masters._masters_hash`) and a reference solution's hash
(:func:`~mastpipe.steps.wavecal.solutions_hash`)."""
