# Install Python packages not in conda
python -m pip install \
  magna==3.4.1

# Patch GUNC
SITE_PACKAGES=$(python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")
patch "$SITE_PACKAGES"/gunc/gunc.py data/gunc.patch

# patch rq
patch "$SITE_PACKAGES"/rq/defaults.py data/rq_defaults.patch
