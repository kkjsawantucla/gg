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

Using ``vib_corrections``
-------------------------

You can pass ``vib_corrections`` to include approximate vibrational free-energy
corrections for adsorbates while constructing phase entries. The value should be
a dictionary mapping species labels (for example ``"OH"`` or ``"H"``) to configured
modifier objects.

Example:

.. code-block:: python

    from gg.modifiers import Remove
    from gg.predefined_sites import FlexibleSites
    from gg.phase_diag import plot_phase_diagram_from_run

    FS2 = FlexibleSites(constraints=True, com=0.75)
    vibOH = Remove(FS2, "OH", weight=0.20)
    vibH = Remove(FS2, "H", weight=0.15)

    fig = plot_phase_diagram_from_run(
        "H2O",
        "H2",
        limits=[[-15, -10], [-10, -5]],
        base_folders="./Mie/",
        mu_path="./Mie/input.yaml",
        file_type=["opt.log", "CONTCAR"],
        vib_corrections={"OH": vibOH, "H": vibH},
        annotate=True,
    )
