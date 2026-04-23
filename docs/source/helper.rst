Helper Utilities
================

Plotting 2D Phase Diagrams
--------------------------

You can generate a 2D phase diagram directly from a completed GCBH run using
``plot_phase_diagram_from_run`` from :mod:`gg.phase_diag`.

This helper scans one or more run folders, collects relaxed structures and energies,
constructs phase entries using the chemical potentials in ``input.yaml``, builds the
2D convex-hull domains, and plots the stable regions.

.. code-block:: python

    from gg.phase_diag import plot_phase_diagram_from_run

    plot_phase_diagram_from_run(
        n1="O",                          # x-axis chemical potential species
        n2="H",                          # y-axis chemical potential species
        limits=[[-3.5, -1.5], [-1.0, 0.5]],
        base_folders=["./"],             # one or more GCBH run directories
        mu_path="./input.yaml",          # must contain chemical_potential
        file_type=["OSZICAR", "CONTCAR"],
        read_from_file=False,             # or path like "entries_O_H.json" to reuse entries
        annotate=True,
        number_labels=True,
        vib_corrections=None,
    )

Typical workflow:

1. Run GCBH and ensure relaxed outputs exist in your run folders.
2. Make sure ``input.yaml`` includes ``chemical_potential`` entries for the species used in ``n1`` and ``n2``.
3. Call ``plot_phase_diagram_from_run(...)`` with the correct folders and file names.
4. The function writes cached entries to ``entries_<n1>_<n2>.json`` and then plots the phase diagram.

Notes:

- ``base_folders`` can include multiple runs to combine structures into one diagram.
- ``file_type`` should match how energies/structures are stored in your calculation outputs.
- If you already generated entries once, set ``read_from_file`` to that JSON file path for faster re-plotting.
