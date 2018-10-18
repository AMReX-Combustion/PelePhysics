  setenv DEVELOPER marc

  setenv DV_DIR ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia
  setenv BLD_ROOT ${DV_DIR}/builds
  setenv BLD_CONFIG ${DV_DIR}/config
  setenv PYTHIA_DIR ${DV_DIR}/pythia

  setenv EXPORT_ROOT ${DV_DIR}/pythia-0.4
  setenv PYTHONPATH ${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver

# Python
  setenv PYTHON_VERSION 2.7
  setenv FUEGO_PYTHON python${PYTHON_VERSION}

  rehash
