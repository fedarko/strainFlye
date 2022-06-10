import os
import shutil
import tempfile
import strainflye.misc_utils as mu

def test_make_output_dir_exists():
    with tempfile.TemporaryDirectory() as tdname:
        assert os.path.exists(tdname)
        child_file_name = os.path.join(tdname, "test_file.txt")
        with open(child_file_name, "w") as f:
            f.write("please don't overwrite me!")

        mu.make_output_dir(tdname)

        # Check that this didn't overwrite anything (the directory itself, or
        # any of its contents)
        assert os.path.exists(tdname)
        assert os.path.exists(child_file_name)
        with open(child_file_name, "r") as f:
            ctext = f.read()
        assert ctext == "please don't overwrite me!"


def test_make_output_dir_notexists_multilevel():
    tdname = os.path.join(tempfile.gettempdir(), "_strainflye_td1", "_strainflye_td2")
    try:
        assert not os.path.exists(tdname)
        mu.make_output_dir(tdname)
        assert os.path.exists(tdname)
    finally:
        # https://stackoverflow.com/a/13118112
        shutil.rmtree(tdname)
