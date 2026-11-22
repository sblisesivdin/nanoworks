Slurm: submit a calculation
===========================

.. literalinclude:: ../../../nanoworks/examples/slurm-profiles/truba-example.json
   :language: json
   :caption: truba-example.json

``time`` limits runtime; ``memory`` requests RAM. ``partition``, ``account``
and ``qos`` select scheduler settings. ``modules`` loads the QE environment;
``job_name`` names the job.
Replace account, partition, QoS, module and executable values with those
of your cluster. They are examples, not Nanoworks requirements.

Run
---

After editing the profile for your cluster, prepare a QE silicon job:

.. code-block:: bash

   cd ~/.nanoworks/examples/Si-2atoms-optical
   dftsolve -E QE -p 4 -i Si-Step1-ground_dos_band.py \
     -g Si_mp-149_primitive_Example.cif -s DOS_integration=tetrahedron \
     --dry-run --cluster-profile ../slurm-profiles/truba-example.json

Submit the generated ``.slurm`` file with ``sbatch`` from its result directory.
The profile supplies scheduler settings; the Python input supplies the physics.

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

Read the submitted job's status and logs, then interpret its numerical
results as for a local calculation. Submission alone does not establish
that the calculation completed. The bundled profile README describes
each field.
