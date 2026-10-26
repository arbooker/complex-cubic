#!/usr/bin/env python3
"""Tests for release-driver fail-closed behavior and manifest checking."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import check_data
import run_certify

RECORD = "-199|1,2,3,1|1|1,0,0,0,1,0,0,0,1|0,1,0|n 1 1 0 0"
VERDICT = "CERTIFIED m = 2 >= |N beta| = 2"


class ReleaseTests(unittest.TestCase):
    def test_threshold_and_exact_outputs(self):
        self.assertEqual(run_certify.parse_output(VERDICT, [RECORD], False)[0]["status"], "certified")
        row = run_certify.parse_output("# comment\nMINIMAL 1 0 -2\n" + VERDICT, [RECORD], True)[0]
        self.assertEqual(row["minimum_coordinates"], [1, 0, -2])

    def test_missing_extra_or_bad_verdict_fails(self):
        for output, exact in [("", False), (VERDICT + "\n" + VERDICT, False),
                              ("REFUSED failed", False), ("REFUTED failed", False),
                              (VERDICT, True), ("CANDIDATE 1 0 0\n" + VERDICT, True),
                              ("MINIMAL 1 0 0\nCERTIFIED lower bound only", True)]:
            with self.subTest(output=output):
                with self.assertRaises(ValueError):
                    run_certify.parse_output(output, [RECORD], exact)

    def test_nonzero_exit_rejected_even_with_valid_output(self):
        result = subprocess.CompletedProcess([], 1, VERDICT, "error")
        with patch.object(run_certify.subprocess, "run", return_value=result):
            with self.assertRaises(RuntimeError):
                run_certify.run_chunk("certify", [RECORD], False)

    def test_manifest_detects_changed_missing_and_extra_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            item = root / "source.txt"
            item.write_text("original")
            check_data.manifest(root, True)
            self.assertEqual(check_data.manifest(root), 1)
            item.write_text("changed")
            with self.assertRaises(ValueError):
                check_data.manifest(root)
            item.unlink()
            with self.assertRaises(ValueError):
                check_data.manifest(root)
            item.write_text("original")
            (root / "extra.txt").write_text("extra")
            with self.assertRaises(ValueError):
                check_data.manifest(root)


if __name__ == "__main__":
    unittest.main()
