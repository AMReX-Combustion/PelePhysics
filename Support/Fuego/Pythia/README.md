# Pythia

## Recommended usage

This usage ensures that all necessary python packages are installed
and configured properly..

1. Install the required software
- [pyenv](https://github.com/pyenv/pyenv#installation)
- [poetry](https://python-poetry.org/docs/#installation)

2. Install the environment
```shell
$ poetry env use -- $(which python)
$ poetry install
```

3. Generate a mechanism with
```shell
$ poetry run bash -c "cd ../Mechanism/Models/air && ./make-mechanism.sh"
```
or regenerate all mechanisms with:
```shell
$ poetry run bash -c "cd ../Mechanism/Models/ && python ./useful_script/script_generate_mechs.py useful_script/list_mech"
```

## Other usage

If you have a python environment where the necessary packages are
installed, you can generate a mechanism with:

```shell
$ # source your python environment
$ cd ../Mechanism/Models/air
$ ./make-mechanism.sh
```
