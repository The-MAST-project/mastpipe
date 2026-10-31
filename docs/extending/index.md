# Extending the pipeline

These pages are for the person who will change the code: a student adding a
grating or an instrument, a developer adding a step or a product. They explain
the architecture and its one strict rule, the data model that the layers
exchange, how to add configuration without code, how the simulator and tests
work, and the conventions the repository enforces.

```{toctree}
:maxdepth: 1

architecture
data-model
adding-a-grating
adding-an-instrument
simulator-and-tests
conventions
```

The [API reference](../api/index.md) is generated from the docstrings, module
by module.

## Development setup

```bash
git clone <repository>
cd mastpipe
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs]"
pytest -m "not slow"                                  # fast feedback
ruff check src tests scripts && ruff format src tests scripts
python scripts/privacy_scan.py
sphinx-build -W -b html docs docs/_build/html         # this documentation
```

The repository's `CLAUDE.md` states the working rules for humans and agents
alike: plan before substantial work, test first, keep `DECISIONS.md` current,
never commit private data or local paths.
