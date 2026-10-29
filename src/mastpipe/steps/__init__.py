"""Step modules: pure functions over arrays, tables and plain parameter values.

A step never reads a file, never sees an :class:`~mastpipe.config.models.EffectiveConfig`
and never logs; the pipeline layer feeds it what it needs and writes what it returns.
"""
