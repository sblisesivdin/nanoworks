# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

__version__ = "27.2.0rc0"


def format_startup_banner(command):
    """Use one concise versioned startup line across calculation commands."""
    return f'Nanoworks {__version__} | {command} started'
