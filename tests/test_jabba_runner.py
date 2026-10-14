import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_jabba import convergence_status, native_options


class JaBbARunnerTests(unittest.TestCase):
    def test_default_options_are_explicit_and_use_gurobi(self):
        args, options = native_options({"memory_gb": 64, "options": {}}, "", "", 8)
        self.assertEqual(options["field"], "ratio")
        self.assertEqual(args[args.index("--gurobi") + 1], "TRUE")
        self.assertEqual(args[args.index("--mem") + 1], "51")
        self.assertNotIn("--purity", args)

    def test_partial_metadata_fails_if_enabled(self):
        with self.assertRaisesRegex(ValueError, "both"):
            native_options({"memory_gb": 64, "options": {},
                            "use_metadata_purity_ploidy": True}, 0.6, "", 8)

    def test_cplex_selection_omits_gurobi_flag(self):
        args, _ = native_options({"memory_gb": 64, "options": {},
                                  "solver": "cplex"}, "", "", 12)
        self.assertNotIn("--gurobi", args)
        self.assertEqual(args[args.index("--cores") + 1], "1")

    def test_subgraph_parallelism_is_explicit(self):
        args, _ = native_options({"memory_gb": 64, "options": {},
                                  "solver": "cplex", "subgraph_cores": 2}, "", "", 12)
        self.assertEqual(args[args.index("--cores") + 1], "2")

    def test_pinned_solver_thread_budget_cannot_oversubscribe(self):
        with self.assertRaisesRegex(ValueError, "Allocate at least"):
            native_options({"memory_gb": 64, "options": {}, "solver": "cplex",
                            "solver_threads": 12, "subgraph_cores": 2}, "", "", 12)

    def test_unknown_solver_fails(self):
        with self.assertRaisesRegex(ValueError, "solver"):
            native_options({"memory_gb": 64, "options": {},
                            "solver": "other"}, "", "", 8)

    def test_unknown_native_option_fails(self):
        with self.assertRaisesRegex(ValueError, "Unknown"):
            native_options({"memory_gb": 64, "options": {"made_up": 1}}, "", "", 8)

    def test_nonconverged_feasible_solution_is_not_success(self):
        result = convergence_status("0.01", 1e-5,
                                    "JaBbA did not converge. Try increasing the tilim.")
        self.assertFalse(result["converged"])
        self.assertEqual(result["status"], "review_required_nonconverged")

    def test_converged_solution_is_success(self):
        result = convergence_status("0", 1e-5, "Done .. job output in: /tmp/output")
        self.assertTrue(result["converged"])

    def test_missing_completion_marker_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "completion marker"):
            convergence_status("0", 1e-5, "unexpected end")


if __name__ == "__main__":
    unittest.main()
