## Contributing

```bash
pip install -e ".[dev,plot]"
pytest              # unit tests, no external dependencies
ruff check .
```

Oracle-comparison tests (`pytest -m oracle`, and the wider set of tests that
compare against the legacy C for extra confidence) need the RUMP C source,
which isn't redistributed here — see [Design and validation](../physics/validation.md).
They skip cleanly when it's absent, so it's not needed for everyday development.
