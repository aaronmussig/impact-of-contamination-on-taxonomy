# Patch GUNC
SITE_PACKAGES=$(python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")
patch "$SITE_PACKAGES"/gunc/gunc.py data/gunc.patch

# patch rq
patch "$SITE_PACKAGES"/rq/defaults.py data/rq_defaults.patch
