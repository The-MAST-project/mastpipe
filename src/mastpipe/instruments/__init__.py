"""Instrument definitions as package data, located and loaded by
:mod:`mastpipe.instruments.registry`.

Each instrument is a folder holding ``config.yaml`` (cameras, step parameters, layouts),
``header_map.yaml`` (raw-header vocabulary), one YAML file per grating or band under
``gratings/`` or ``bands/``, fiber layouts under ``layouts/``, line tables under ``lines/``
and anchor pairs under ``anchors/``. Nothing in this package is code that branches on an
instrument or a grating name."""
