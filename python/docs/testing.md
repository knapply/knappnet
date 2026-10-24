# Testing

## The backend's own tests

```sh
python -m pip install "./python[test]"
python -m pytest python/tests
```

`python/tests/test_against_networkx.py` calls every implemented function
with `backend="knappnet"` and with `backend="networkx"` on several hundred
graph and argument combinations and compares the results. It fails if any
implemented function never ran on the engine.

## NetworkX's test suite

NetworkX can run its whole test suite through a backend: each dispatchable
call converts its graphs to the backend and calls the backend's
implementation. With fallback enabled, functions knappnet does not run go
to NetworkX:

```sh
python -m pip install pytest numpy scipy
NETWORKX_TEST_BACKEND=knappnet NETWORKX_FALLBACK_TO_NX=True \
    python -m pytest --pyargs networkx -p no:cacheprovider
```

Tests that check behaviour knappnet deliberately does not share (see
[differences](differences.md)) are marked as expected failures by the
backend's `on_start_tests` hook (`knappnet/_testing.py`).
