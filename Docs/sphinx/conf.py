extensions = [ 'sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'PelePhysics'
copyright = u'AMReX Copyright (c) 2022, The Regents of the University of California, through Lawrence Berkeley National Laboratory and the Alliance for Sustainable Energy, LLC., through National Renewable Energy Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.'
author = u'J.B. Bell, M.S. Day, E. Motheau, D. Graves, M. Henry de Frahan, R.W. Grout'
version = u'2022.10'
release = u'2022.10'
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}
html_theme = 'sphinx_rtd_theme'
htmlhelp_basename = 'PelePhysicsdoc'
latex_elements = {
}
latex_documents = [
    (master_doc, 'PelePhysics.tex', u'PelePhysics Documentation',
     author, 'manual'),
]
man_pages = [
    (master_doc, 'pelec', u'PelePhysics Documentation',
     [author], 1)
]
texinfo_documents = [
    (master_doc, 'PelePhysics', u'PelePhysics Documentation',
     author, 'PelePhysics', 'One line description of project.',
     'Miscellaneous'),
]
