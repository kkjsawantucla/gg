GCBH
====

.. contents::
    :local:

Calculator
----------
The ``Gcbh`` calculator extends ASE Dynamics to run a grand-canonical basin-hopping workflow for
surface/cluster structure search. At each step it applies a structure modifier, performs local
optimization, evaluates grand-canonical free energy using user-provided chemical potentials, and
accepts/rejects moves with a Metropolis criterion.

The run writes standard outputs such as ``gcbh.log`` (run history), ``gcbh.traj`` (trial
structures), ``local_minima.traj`` (accepted minima), and ``current_status.pkl`` (restart state),
with optimization artifacts in ``opt_folder/``.

Most runtime functionality is controlled through ``input.yaml`` (loaded via ``config_file``),
including temperature schedule, stopping rules, optimization behavior, uniqueness checks, and
accept/reject safety checks.

For examples, see the `gg/examples/ <https://github.com/kkjsawantucla/gg/blob/main/examples/>`_
folder or the detailed explanation in :ref:`Examples`.

Attaching an ASE calculator (MACE, FAIR-Chem, NequIP)
-----------------------------------------------------

``Gcbh`` requires ``atoms.calc`` to be set before initialization. Any ASE-compatible calculator
can be used. Below are common ML calculator setups.

MACE
~~~~

.. code-block:: python

    from ase.io import read
    from mace.calculators import mace_mp

    atoms = read("POSCAR")
    calc = mace_mp(model="medium-mpa-0", device="cuda")  # or device="cpu"
    atoms.calc = calc

Meta FAIR-Chem (Open Catalyst / fairchem)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from ase.io import read
    from fairchem.core.common.relaxation.ase_utils import OCPCalculator

    atoms = read("POSCAR")
    calc = OCPCalculator(
        model_name="eqv2_31M_omol",  # replace with your model name
        checkpoint_path="path/to/checkpoint.pt",  # optional if model is cached
        cpu=False,
    )
    atoms.calc = calc

NequIP
~~~~~~

.. code-block:: python

    from ase.io import read
    from nequip.ase import NequIPCalculator

    atoms = read("POSCAR")
    calc = NequIPCalculator.from_deployed_model("path/to/nequip_model.pth")
    atoms.calc = calc

After attaching the calculator, create ``Gcbh(atoms, config_file="input.yaml")`` and add
modifiers before calling ``run(...)``.

.. autoclass:: gg.gcbh.Gcbh
    :members: run
    :undoc-members:
    :show-inheritance:

.. autoclass:: gg.gcbh.GcbhFlexOpt
    :members: run
    :undoc-members:
    :show-inheritance:

Attaching modifiers to GCBH
---------------------------

``Gcbh`` needs at least one structure modifier before ``run(...)`` can sample new structures.
Attach modifiers with ``add_modifier(instance, name)`` where:

* ``instance`` is a modifier object (for example ``Add``, ``Remove``, ``Replace``, ``Swap``, or ``ModifierAdder``).
* ``name`` is a unique label used in logs and for adaptive modifier weighting.

If a duplicate ``name`` is provided, ``Gcbh`` raises an error.

.. code-block:: python

    from gg.gcbh import Gcbh
    from gg.modifiers import Add, Remove

    G = Gcbh(atoms, config_file="input.yaml")

    add_o = Add(surface_sites, "O", surf_coord=[2, 3], surf_sym=["Cu"])
    rem_o = Remove(cluster_sites, "O")

    G.add_modifier(add_o, "Add O")
    G.add_modifier(rem_o, "Remove O")

    # Optional: remove detached gas fragments during acceptance
    G.add_delete_gas(gas_species=["H2"])

    G.run(steps=1000)

input.yaml configuration reference
----------------------------------

The ``config_file`` argument in ``Gcbh(...)`` loads key/value pairs from ``input.yaml`` and
updates the internal ``self.c`` dictionary. Keys below are the primary runtime controls and
their defaults in :class:`gg.gcbh.Gcbh`.

.. code-block:: yaml

    temp: 1500
    max_temp: null
    min_temp: null
    stop_steps: 40
    stop_opt: 500
    vasp_opt: false
    chemical_potential: null
    max_history: 25
    max_bond: 2
    max_bond_ratio: 0
    check_graphs: true
    area: false
    fmax: 0.05
    opt_traj: false
    vib_correction: false
    initialize: true
    detect_gas: null
    graph_method: fullgraph
    check_contact_error: 0.5

Parameter details:

* ``temp`` (float): Metropolis temperature (K-like scale used with ``kB`` in acceptance probability).
* ``max_temp`` (float or null): Upper bound for adaptive temperature. If ``null``, set to ``1.5 * temp`` at startup.
* ``min_temp`` (float or null): Lower bound for adaptive temperature. If ``null``, set to ``temp / 1.5`` at startup.
* ``stop_steps`` (int): Stop basin-hopping after this many **consecutive** non-improving accepted/rejected steps.
* ``stop_opt`` (int): Maximum geometry-optimization steps per trial move.
* ``vasp_opt`` (bool): If ``true``, skip ASE optimizer and only evaluate potential energy/contact checks (useful with external workflows).
* ``chemical_potential`` (dict): Required mapping like ``{"Cu": -4.1, "O": -7.5}`` used to compute grand-canonical free energy ``F = E - sum(mu_i * n_i)``.
* ``max_history`` (int): Length of acceptance history window used to adapt temperature.
* ``max_bond`` (float): Bond-length scaling factor used in graph construction/removal-site logic.
* ``max_bond_ratio`` (float): Additional bond tolerance used by graph/modifier connectivity checks.
* ``check_graphs`` (bool): If ``true``, reject duplicate structures based on graph uniqueness checks.
* ``area`` (bool): If ``true``, report/track energy and free energy normalized by surface area after initialization.
* ``fmax`` (float): Force convergence criterion for ASE optimizer.
* ``opt_traj`` (bool or int): If truthy, write ``opt.traj`` during optimization; integer values are used as write interval.
* ``vib_correction`` (bool): Enable optional vibrational correction hooks added via ``add_vib_correction(...)``.
* ``initialize`` (bool): If ``false``, skip initial optimization and start from current structure with placeholder thermodynamics.
* ``detect_gas`` (null/list): Reserved key in current implementation; gas removal is typically controlled via ``add_delete_gas(...)``.
* ``graph_method`` (str): Graph-comparison strategy passed to uniqueness checks. Use ``"fullgraph"`` for whole-structure matching, or provide an atom symbol (for example ``"O"`` or ``"Cu"``) to build centered subgraphs from that species when checking uniqueness.
* ``check_contact_error`` (float): Contact threshold passed to ``check_contact(...)``; structures with touching atoms are rejected.
