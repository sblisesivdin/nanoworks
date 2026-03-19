Contributing
============

We welcome contributions and hope this guide will make the `Nanoworks` code repository easier to understand. It is important to mention that the `Nanoworks` software development is run voluntarily and therefore we need to build a community that can support user questions, attract new users, maintain documentation, write tutorials, and develops new features to make this software a useful tool for all users.

Being Respectful
----------------
Please show empathy and kindness towards other people, software, and communities working diligently to develop other tools.

Please do not speak negatively about other people or their work in pull requests and issues.

Cloning the Source Repository
-----------------------------

Before cloning the source repository to your computer, please visit the [installation page](https://sblisesivdin.github.io/nanoworks/installation/) of `Nanoworks` to install ASE, GPAW, and all other needed packages to your computer. Then, you can clone the source of the `Nanoworks` from the related repository:
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

New Release
-----------

In each new release, the following steps must be completed.
- Update `__version__` variable in `‎nanoworks/__init__.py` as YY.m.x. Here, YY is the year, m (or mm) is the month, and x is the step number (first step is 0) of releases for that month. There is no minor or major release.
- Update `‎RELEASE_NOTES.md` and `‎pyproject.toml` as done in this [commit](https://github.com/sblisesivdin/nanoworks/commit/0593543fe472815dcd0736175b7e5249eeffa74c#diff-68d24fa81558aae3d8c59e2aa57a4fa719ea3b04d7fa14beff45f16f00858f50).
- Update `release` and `version` in `‎‎docs/source/conf.py` as done in this [commit](https://github.com/sblisesivdin/nanoworks/commit/1f460ceca08062b16590d2010ac4c3f1194f6a36).
- On the releases page, click `Draft a new release` and create a new tag with the version number.
- Give a general title, give some highlights information, write Release Notes (copy/paste from RELEASE_NOTES.md page, it can be automatized in future).
- Select `Create a discussion for this release` and then finish the release.
