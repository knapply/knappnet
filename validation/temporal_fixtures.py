"""Contact times for the temporal-network validation.

Temporal measures are validated on the ordinary fixture corpus
(validation/fixtures/), with every edge read as a contact whose start and
end times follow deterministically from its position in the edge list.
Times are small integers, so they are exact in floating point and many
contacts share a start time, which exercises the rules for simultaneous
contacts; about 40% of contacts are instantaneous.

Both the engine runners (validation/bindings/temporal.py) and the oracles
take times from here, so they always see the same contacts.
"""

from __future__ import annotations

DURATIONS = (0, 1, 0, 2, 3, 0, 1)


def period(m: int) -> int:
    """Number of distinct start times: about three contacts per time."""
    return max(3, m // 3 + 2)


def contact_times(fx, instantaneous: bool = False) -> tuple[list[float], list[float]]:
    """Start and end time of every edge of a fixture, in edge order."""
    m = len(fx.src)
    p = period(m)
    start = [float((e * 37 + 11) % p) for e in range(m)]
    if instantaneous:
        return start, list(start)
    end = [s + DURATIONS[e % len(DURATIONS)] for e, s in enumerate(start)]
    return start, end
