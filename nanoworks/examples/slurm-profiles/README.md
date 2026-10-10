# Slurm cluster profile example

`truba-example.json` is a template for QE dry-run/Slurm execution. Replace
`YOUR_PROJECT`, partition, QoS and module names with the values provided by
your cluster. These example values are not verified TRUBA site settings.

Install the bundled examples with:

```bash
nanoworks --install-examples
```

Edit the profile, then select it directly from your calculation directory:

```bash
dftsolve --dry-run -p 48 \
  --cluster-profile ~/.nanoworks/examples/slurm-profiles/truba-example.json \
  -i input.py -g geometry.cif
```

The input must select QE (`Engine = 'QE'` or `-E QE`). The command prepares
input files and a Slurm script; it does not submit a job. Review the generated
script before submitting it with `sbatch`. For elasticity, also load the
cluster's thermo_pw 2.1.0 module if it is separate from QE.

Alternatively, copy the edited profile to
`~/.config/nanoworks/clusters/truba.json` and use `--cluster-profile truba`.
Explicit `--slurm-*` options override profile values.

`--install-examples` preserves an existing examples directory. If you already
installed examples before this profile was bundled, copy this folder from the
source checkout into your existing examples directory.
