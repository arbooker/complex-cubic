#!/usr/bin/env python3
"""Regression tests for verification bookkeeping; no proof binaries are run.

    python3 -B src/test_verification.py
"""
import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import cmp_nfdata
import run_euclid


FIELD = "-23|1,0,-1,1|1|1,0,0,0,1,0,0,-1,1|0,1,0|-23\n"
SECOND = "-31|-1,1,0,1|1|1,0,0,1,0,1,0,1,0|0,1,0|-31\n"


class VerificationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="verification-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def compare(self, first, second):
        a, b = self.root / "a.dat", self.root / "b.dat"
        a.write_text(first)
        b.write_text(second)
        with redirect_stdout(io.StringIO()):
            return cmp_nfdata.main(str(a), str(b))

    def test_matching_bases_pass(self):
        self.assertEqual(self.compare(FIELD + SECOND, FIELD + SECOND), 0)

    def test_missing_record_in_either_file_fails(self):
        self.assertNotEqual(self.compare(FIELD + SECOND, FIELD), 0)
        self.assertNotEqual(self.compare(FIELD, FIELD + SECOND), 0)

    def test_empty_and_disjoint_comparisons_fail(self):
        for a, b in [("", ""), (FIELD, ""), ("", FIELD), (FIELD, SECOND)]:
            with self.subTest(a=a, b=b):
                self.assertNotEqual(self.compare(a, b), 0)

    def test_fixed_list_is_complete_and_unambiguous(self):
        discs = run_euclid.EUCLIDEAN
        self.assertEqual(len(discs), 58)
        self.assertEqual(len(set(discs)), 58)
        self.assertEqual(set(run_euclid.field_records(discs)), set(discs))

    def test_missing_or_ambiguous_field_records_fail(self):
        nf = self.root / "nf.dat"
        with patch.object(run_euclid, "NFDATA", str(nf)):
            nf.write_text(FIELD)
            with self.assertRaisesRegex(ValueError, "missing"):
                run_euclid.field_records([-23, -31])
            nf.write_text(FIELD + FIELD)
            with self.assertRaisesRegex(ValueError, "ambiguous"):
                run_euclid.field_records([-23])

    def invoke_driver(self, fail=False):
        # Use the actual default output path, not --out, to exercise the bug
        # where the output used to be read back as the next run's work list.
        nf = self.root / "nf.dat"
        nf.write_text(FIELD + SECOND)
        seen = []

        def fake_run(binary, path, depth, gbox):
            disc = int(Path(path).read_text().split("|", 1)[0])
            seen.append(disc)
            # Two ordering calls, then a completed result, then interruption.
            if fail and len(seen) == 4:
                raise KeyboardInterrupt()
            return 0.01, "cubes=10 deepest=2\nPROVED every point covered\n"

        argv = ["run_euclid.py", "--jobs", "1", "--tmp", str(self.root)]
        with patch.object(run_euclid, "ROOT", str(self.root)), \
             patch.object(run_euclid, "NFDATA", str(nf)), \
             patch.object(run_euclid, "EUCLIDEAN", [-23, -31]), \
             patch.object(run_euclid, "run", side_effect=fake_run), \
             patch("sys.argv", argv), redirect_stdout(io.StringIO()):
            result = run_euclid.main()
        self.assertFalse(list(self.root.glob("euclid-*")))
        self.assertEqual(seen, [-23, -31, -23, -31])
        return result

    def output_path(self):
        results = self.root / "work"
        results.mkdir(exist_ok=True)
        return results / "euclidean.jsonl"

    def test_partial_previous_output_does_not_shrink_work_list(self):
        out = self.output_path()
        out.write_text(json.dumps({"disc": -23, "status": "proved"}) + "\n")
        self.assertEqual(self.invoke_driver(), 0)
        records = [json.loads(s) for s in out.read_text().splitlines()]
        self.assertEqual([r["disc"] for r in records], [-23, -31])

    def test_no_previous_output_is_required(self):
        out = self.output_path()
        self.assertEqual(self.invoke_driver(), 0)
        self.assertEqual(len(out.read_text().splitlines()), 2)

    def test_interruption_preserves_old_results(self):
        out = self.output_path()
        original = b"previous complete results\n"
        out.write_bytes(original)
        with self.assertRaises(KeyboardInterrupt):
            self.invoke_driver(fail=True)
        self.assertEqual(out.read_bytes(), original)
        self.assertEqual(list(out.parent.iterdir()), [out])
        self.assertFalse(list(self.root.glob("euclid-*")))

    def test_failed_replace_preserves_old_results(self):
        out = self.output_path()
        out.write_text("old\n")
        with patch.object(os, "replace", side_effect=OSError("replace failed")):
            with self.assertRaises(OSError):
                with run_euclid.atomic_results(out) as pending:
                    pending.write("new\n")
        self.assertEqual(out.read_text(), "old\n")
        self.assertEqual(list(out.parent.iterdir()), [out])


if __name__ == "__main__":
    unittest.main()
