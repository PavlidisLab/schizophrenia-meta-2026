# provenance/

Tables and scripts cited in the hidden source comments of `manuscript_ph.md`, `supplement.md` and `Figures_ph.md`,
copied from the analysis repository and deliveries. In the copies, home-directory and working-folder path prefixes,
host names and credential lookups were removed; the content is otherwise unchanged. One folder per place the files are
cited: `main_text`, `extended_methods_NN_*`, `extended_figure_NN`, `supplementary_note_NN_*`, `supplementary_table_N`
and `supplementary_data_1`, with scripts in a `scripts/` subfolder where a folder has them. Folder numbers need not match
the current section numbering. A file cited in several places sits under
the first. `MANIFEST.tsv` gives the sha256 of each file; files too large for the repository are marked `to_deposit`.
Files cited but not distributed are marked in the comments `(not in the repository; see DATA_NEEDED.tsv)` or
`(script not distributed)`; `tools/check_source_comments.py` checks that every cited path is in the repository or
marked.
