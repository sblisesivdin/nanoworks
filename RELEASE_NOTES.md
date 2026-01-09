## Release Notes

### Development Version

- Initializing context and tools for the session.
- Removed gg.py and gui_files/ directory.
- Renaming project and its files as: gpaw-tools -> Nanoworks, gpawsolve.py -> dftsolve.py, and asapsolve -> mdsolve.py
- Replace the gpaw-tools setup with the Nanoworks setup
- Implement an ML solver script: mlsolve.py
- Remove GW calculations and related parameters from dftsolve.py
- Simplify GPW file writing logic (Always mode="all")
- A major refactorization of global variable usage to dataclass! This change ensures type safety, centralized configuration, reduced global namespace pollution, and improved code organization with total backward compatibility.
- In addition to the refactorization of global variable usage, more than 30 security warnings have been fixed.
- An easy example for machine learning capabilities is added as examples/Graphene-ML
- Behavior of output folder creation and writing in that directory is now same for dftsolve.py, mdsolve.py and mlsolve.py
- Drawing to screen -d argument is removed.
- pip install nanoworks works.
- Website is moved to https://nanoworks.readthedocs.io
- dftsolve.py: work as system-wide command "dftsolve"
- mdsolve.py: work as system-wide command "mdsolve"
- mlsolve.py: work as system-wide command "mlsolve"
- New command: "nanoworks" command is added. Just giving important information for now.
