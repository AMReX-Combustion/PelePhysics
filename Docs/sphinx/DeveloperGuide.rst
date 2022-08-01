.. highlight:: rst


Developer Guidelines
====================

The following developer guidelines apply to all code that is to be committed to the `Pele` suite of codes. 

Before adding files for a commit to be pushed to a branch of `PelePhysics`, first run the following formatting tools depending on the type of code:


C++ Code
--------
All C++ code should be processed by the Clang formatter prior to being added for commit.

`clang-format -i .` Will apply all of the correct formatting conventions to make all code consistent.


Python
------

The tools necessary to format the `Python` code are maintained through `poetry`.

1) ``isort``: `poetry run isort .` This will perform proper sorting of the installed `Python` libraries.
2) ``black``: `poetry run black --preview .` This will perform proper formatting of all `Python` files to be consistent with all current files.
3) ``flake8``: `poetry run flake8 .` This will run diagnostics on all the `Python` files and will list a series of issues that need to be addressed to adhere to current `Python` best practices


