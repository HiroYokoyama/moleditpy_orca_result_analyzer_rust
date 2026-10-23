# Source and license notices

The shared Python UI, analysis modules, parser compatibility mixins and tests
were synchronized from HiroYokoyama/moleditpy_orca_result_analyzer_plugin at
commit d7f79374829da98c2f4c84dc6811504d0f5fc164 (plugin 5.0.3).
Copyright: Hiromichi Yokoyama / HiroYokoyama and the existing project notices.
That upstream source is available under LGPL-3.0-or-later; retain the included
licenses/LGPL-3.0.txt and GPL text in LICENSE with these copies.
The existing GPLv3 license of the Rust edition is retained; this synchronization
does not change third-party or native dependency licenses.

The native Rust parser and MO evaluator remain local. The Python parser mixins
supply shared helper APIs and NBO/spin compatibility parsing; the main numerical
MO grid evaluation and bulk output parsing run in Rust.

Rust dependencies (including PyO3, regex, once_cell and Rayon) retain their own
licenses. Python, Qt, RDKit, NumPy, Matplotlib and optional runtime dependencies
also retain their own terms. Test data provenance is recorded in
tests/sample_outputs/NOTICE.
