import os
import sys
import pathlib
import contextlib
from shutil import copytree, rmtree
import ceptr.ceptr as ceptr


mechanisms = ['grimech12', 'Davis', 'heptane_3sp', 'drm19', 'heptane_fc', 'LiDryer', 'SootReaction',
            'grimech30-noArN', 'air', 'grimech30', 'alzeta', 'propane_fc', 'ethylene_af', 'Kolla',
            'heptane_lu_88sk', 'LuDME', 'nitrogens', 'sCO2', 'dodmethair_4sp', 'JL4', 'chem-H',
            'BurkeDryer', 'FFCM1_Red', 'decane_3sp', 'dodecane_lu', 'dodecane_wang', 'chem-CH4-2step']


class DummyFile(object):
    def write(self, x): pass


@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = DummyFile()
    yield
    sys.stdout = save_stdout


def get_raw_file(filepath):
    f = open(filepath, 'r')
    f_content = f.read()
    f_content = f_content.replace(' ', '')
    f_content = f_content.replace('\n', '')
    return f_content


def mechanism_path(mname):
    """Determine mechanism path."""
    this_file_dir = pathlib.Path(__file__).parent.resolve()
    model_path = "Mechanism/Models"
    return this_file_dir.parents[1] / model_path / mname


def test_regression_ceptr():
    """ Test that a new version produces the same mechanism file as before """

    # Loop over all the possible mechanism files
    for i, mech in enumerate(mechanisms):
        # Check that it is a convertible mechanism
        mech_path = mechanism_path(mech)
        fname = mech_path / "mechanism.yaml"

        # Clean up from previous runs
        if 'temp_mech_{i}' in os.listdir():
            rmtree('temp_mech_{i}')

        # Copy the mechanism file to the current directory
        copytree(f'{mech_path}', f'temp_mech_{i}')
        
        # Rename the mechanism files in the mechanism
        os.rename(f'temp_mech_{i}/mechanism.H', f'temp_mech_{i}/mechanism_old.H')
        os.rename(f'temp_mech_{i}/mechanism.cpp', f'temp_mech_{i}/mechanism_old.cpp')
    
        # Run the convertion utility without output
        with nostdout():
            ceptr.convert(f'temp_mech_{i}/mechanism.yaml', None, False)
    
        # Compare new with old
        result_header = get_raw_file(f'temp_mech_{i}/mechanism.H') == get_raw_file(f'temp_mech_{i}/mechanism_old.H')
        result_cpp = get_raw_file(f'temp_mech_{i}/mechanism.cpp') == get_raw_file(f'temp_mech_{i}/mechanism_old.cpp')

        if not result_header:
            rmtree(f'temp_mech_{i}')
            raise(AssertionError(f'Regression test failed for the header file for mechanism {mech})'))

        elif not result_cpp:
            rmtree(f'temp_mech_{i}')
            raise(AssertionError(f'Regression test failed for the cpp file for mechanism {mech})'))
        else:
            rmtree(f'temp_mech_{i}')

                                 
            
