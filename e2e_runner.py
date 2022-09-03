#!/usr/bin/env python3.8
"""
E2E/Snapshot testing tool

Test frame3dd binary against all examples
Checks output, files count and content with reference files
"""

import difflib
import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path


EXECUTABLE_PATH = "build/frame3dd"
WORKDIR_PATH = ".e2e_workdir"
WORKDIR = Path(WORKDIR_PATH)


@contextmanager
def use_chdir(path):
    old_path = os.getcwd()
    try:
        yield os.chdir(path)
    finally:
        os.chdir(old_path)


def check_executable():
    """Check if the code is built"""
    try:
        os.stat(EXECUTABLE_PATH)
    except FileNotFoundError:
        print(
            f"Executable not found at {EXECUTABLE_PATH}. Did you forget to build the code?",
            file=sys.stderr,
        )
        sys.exit(1)


def set_up():
    """Pre-test routines"""
    try:
        os.mkdir(WORKDIR)
    except FileExistsError:
        print(
            f"Warning! Workdir was not removed from previous run. Cleaning it up now",
            file=sys.stderr,
        )
        tear_down()
        os.mkdir(WORKDIR)

    os.symlink("../" + EXECUTABLE_PATH, WORKDIR / "frame3dd")


def tear_down():
    """Post-test routines"""
    shutil.rmtree(WORKDIR)


def test_case(example_name):
    examples = Path("../examples")
    process = subprocess.run(
        ["./frame3dd", "-i", examples / f"{example_name}.3dd", "-o", "out"]
    )
    assert process.returncode == 0, f"Exit code is non-zero! {process.returncode}"
    with open(examples / f"{example_name}.out") as fd:
        original_lines = fd.readlines()
    with open("out") as fd:
        produced_lines = fd.readlines()

    diff = difflib.unified_diff(original_lines, produced_lines, n=0)

    try:
        # Skip junk
        next(diff)
        next(diff)
        while True:
            header = next(diff)
            orig = next(diff)
            actual = next(diff)
            # Ignore data changes
            if header == "@@ -9 +9 @@\n":
                continue
            # FIXME Actually check RMS error
            if orig.startswith("-R M S"):
                continue
            # FIXME Actually check Total Mass deviation
            if orig.startswith("-  Total Mass"):
                continue
            # FIXME Actually check Mode frequency deviation
            if orig.startswith("-  MODE"):
                continue
            raise AssertionError("Regression error:\n" + header + orig + actual)
    except StopIteration:
        pass


check_executable()
for example_name in ("exA", "exB", "exE"):
    set_up()
    with use_chdir(WORKDIR):
        test_case(example_name)
    tear_down()
print("All done")
