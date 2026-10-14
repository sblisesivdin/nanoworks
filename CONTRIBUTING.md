Contributing
============

We welcome contributions and hope this guide will make the `Nanoworks` code repository easier to understand. It is important to mention that the `Nanoworks` software development is run voluntarily and therefore we need to build a community that can support user questions, attract new users, maintain documentation, write tutorials, and develops new features to make this software a useful tool for all users.

Current Development Scope
-------------------------

As of October 11, 2026, development is in feature freeze. Current work focuses
on bug fixes and clearer installation, existing examples and tutorials.
New solver capabilities, input keywords and command-line options are deferred.
Feature ideas may still be recorded for a future development cycle.

Describe how a bug was reproduced and how the correction was checked. Keep
static checks separate from actual engine execution: documentation checks do
not validate GPAW/QE calculations, MPI behavior or numerical accuracy. Before
a release candidate, validate representative calculations on a machine with
the required engines and data, including serial and MPI runs.

Being Respectful
----------------
Please show empathy and kindness towards other people, software, and communities working diligently to develop other tools.

Please do not speak negatively about other people or their work in pull requests and issues.

Cloning the Source Repository
-----------------------------

Use the development route on the [installation page](https://nanoworks.readthedocs.io/en/latest/installation.html) and install dependencies for your selected engine. A QE-only environment does not require GPAW. Clone the source repository:
[Main Github repository](https://github.com/sblisesivdin/nanoworks) by:

    git clone https://github.com/sblisesivdin/nanoworks.git
    cd nanoworks

Questions
---------

We do not have any email list, IRC channel, Slack room, etc. 
* For general questions about the project and all other things, we will use the [Discussions](https://github.com/sblisesivdin/nanoworks/discussions) page of the GitHub repository. 
* For more technical problems, you can create an issue on the [Issues](https://github.com/sblisesivdin/nanoworks/issues) page of the GitHub repository. Posting to the issues page allows community members with the required expertise to answer your question, and the information obtained remains available to other users on the issues page for future usage.

Reporting Bugs and Feature Requests
-----------------------------------

If you encounter any bugs, crashes, or quirks while using the code, please report it on the [Issues](https://github.com/sblisesivdin/nanoworks/issues) page with an appropriate tag so that the developers can take care of it immediately. When reporting an issue, please be overly descriptive so we can reproduce it. Provide trackbacks, screenshots, and sample files to help us resolve the issue. Please create an issue with the "Bug report" template for reporting a bug.

Please do not hesitate to submit ideas for improvements to the `Nanoworks` software. To suggest an improvement, please create an [Issue](https://github.com/sblisesivdin/nanoworks/issues) with the "Feature Request" template. Please use descriptive and extensive information (links, videos, possible screenshots, etc.) to help the developers implement this functionality.

Contributing New Code
---------------------
If you have an idea to solve a bug or implement a new feature, please first create an issue as a bug report/feature request. We can then use that issue as a discussion thread to resolve the contribution implementation.

Licensing
---------

All code is licensed under the MIT License. If you didn't write the code yourself, it's your responsibility to ensure that the existing license is compatible with and included with the contributed files.

The authoritative license text is `LICENSE.md` in the repository root.
Comment-capable source files, tests, executable examples, build scripts and
configuration files carry the following header using the file's comment syntax:

```text
SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
SPDX-License-Identifier: MIT
See LICENSE.md in the project root for license terms.
```

Keep a script's shebang and Python encoding declaration in their required
positions. Preserve existing third-party copyright notices, record included
third-party code in `THIRD_PARTY_NOTICES.md`, and include those notices in
redistributions. Do not insert comments into formats that prohibit them,
such as JSON; these files remain covered by the root license unless otherwise
specified. External dependencies and engines are covered by their own licenses.

Run `python scripts/check_license_headers.py` before submitting new source
files. CI runs the same check without installing or running DFT engines.

New Release
-----------

The next release target is `27.2.0`; the current unreleased source uses `27.2.0rc0`.

In each new release, the following steps must be completed.
- Update `__version__` variable in `‎nanoworks/__init__.py` as YY.m.x. Here, YY is the year, m (or mm) is the month, and x is the step number (first step is 0) of releases for that month. There is no minor or major release.
- Update `RELEASE_NOTES.md` with the final release version and date. Package metadata in `pyproject.toml` reads the version from `nanoworks.__version__`.
- Sphinx `release` and `version` in `docs/source/conf.py` derive from the same package version; verify them instead of maintaining a separate version string.
- On the releases page, click `Draft a new release` and create a new tag with the version number.
- Give a general title, give some highlights information, write Release Notes (copy/paste from RELEASE_NOTES.md page, it can be automatized in future).
- Select `Create a discussion for this release` and then finish the release.
