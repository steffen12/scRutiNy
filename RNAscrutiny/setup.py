from setuptools import setup

setup(
    name='RNAscrutiny',
    version='1.1',
    description='Tool for simulating and analyzing scRNA-seq data.',
    keywords=
    'RNA single cell simulation inference generative genomics transcriptomics',
    url='http://lbm.niddk.nih.gov/mckennajp/scRutiNy',
    author='Steffen K Cromwell, Joseph P McKena, Vipul Periwal',
    author_email='joepatmckenna@gmail.com',
    license='MIT',
    packages=['RNAscrutiny'],
    install_requires=[
        'numpy', 'matplotlib', 'sklearn', 'rpy2', 'seaborn', 'networkx'
    ],
    include_package_data=True,
    zip_safe=False)
