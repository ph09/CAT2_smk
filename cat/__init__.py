"""
CAT pipeline Python package shim.

This repository historically executed modules via `python cat/<script>.py`.
Several scripts also import sibling modules using `from cat.<module> import ...`.
Adding this file makes the local `cat/` directory importable as a Python package.
"""

