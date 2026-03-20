# MLIP: Graphene with Vacancy Example

This example demonstrates how to use `mlsolve` to perform a geometry optimization on a 8x8 graphene structure (128 atoms) and 8x8 graphene structure with a vacancy (127 atoms) using the SevenNet model. You can change model with MACE (mace) and CHGNet (chgnet) in input files.

## Running the Example

To run this example, execute the following command:

    mlsolve -i ml_input-pristine.py -g graphene8x8withvacancy.cif

**NOTE:** This is the example done in Nanoworks article.
