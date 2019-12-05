from setuptools import setup
setup(
    name = "SCGid",
    version = "1.0.0",
    author = "Kevin Amses",
    author_email = "amsesk@umich.edu",
    description = (
        "A consensus approach to contig filtering and genome prediction from single-cell sequencing libraries"
        ),
    license = "GNU General Public License v3.0",
    keywords = "",
    url = "http://www.github.com/amsesk/SCGid",
    packages=["scgid","tests"],
    package_data={
        '':['*.ini', '*.yaml']
    },
    #py_modules=['scripts.scgid', 'scripts.gct', 'scripts.codons', 'scripts.kmers', 'logging'],
    scripts=['bin/scgid'],
    install_requires=[
        'pyyaml',
        'numpy',
        'pandas',
        'ete3'
    ]
)