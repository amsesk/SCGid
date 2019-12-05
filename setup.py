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
    packages=["scgid","scgid.tests"],
    package_data={
        '':['*.ini', '*.yaml']
    },
    #py_modules=['scgid.scgid', 'scgid.gct', 'scgid.codons', 'scgid.kmers', 'logging'],
    scripts=['bin/scgid'],
    install_requires=[
        'numpy>=1.15.0',
        'pandas>=0.23.4',
        'ete3>=3.1.1',
        'plotly',
        'pyyaml'
    ]
)

# Import config script after above
