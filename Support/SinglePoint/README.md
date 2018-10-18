# SinglePoint
* A helper app to enable coupling to PelePhysic modules outside of the standard build system

## Usage

This tool will generate a standalone shared library, libSP.so, to link to applications for evaluating
transport coefficients pointwise.  The library is built in the Exec/Build folder, after setting paths
to the Eos (optionally, Chemistry_Model), Reactions and Transport folders.  A main.cpp is provided
to test the functionality and demonstrate usage.

See the GNUmakefile in the Build folder for more information.

If you run into difficulties, contact Marc Day at MSDay@lbl.gov, or file a GitLab "issue" for support.

