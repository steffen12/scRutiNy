from setuptools import setup

setup(
    name='RNAscrutiny',
    version='0.4',
    description='Generative models from scRNA-seq data.',
    keywords='generative single cell RNA sequencing transcriptomics genomics',
    url='http://www.google.com',
    author='Steffen K Cromwell, Joseph P McKena, Vipul Periwal',
    author_email='joepatmckenna@gmail.com',
    license='MIT',
    packages=['RNAscrutiny'],
    install_requires=[
        'numpy', 'matplotlib', 'sklearn'
    ],
    zip_safe=False)
