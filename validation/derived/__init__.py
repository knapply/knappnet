"""Data derived from the fixture corpus for areas that need more than an
edge list: contact times for temporal networks (temporal.py).

Each module states its rules in its docstring. The engine runner
(validation/bindings/<area>.py), the Python oracles, and the R oracles all
derive the same data from the same rules, so every oracle sees exactly the
same input. Areas whose derived data is simpler (node attributes from the
node index for mixing and brokerage) state the rule in the comment
at the top of their manifest file instead.
"""
